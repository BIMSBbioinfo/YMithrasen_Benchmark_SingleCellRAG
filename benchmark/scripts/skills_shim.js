/**
 * skills_shim.js
 *
 * Sibling of baseline_shim.js, built to SIMULATE the shape of how Claude
 * Skills' progressive disclosure works, as closely as a bare
 * chat-completion API call reasonably allows — not a claim of equivalence
 * to Claude's actual internal skill-loading implementation, which this
 * project has no visibility into. Calls various model providers and
 * returns extracted R code.
 * Supports: gpt-5.5, gpt-5.4-mini, gemini-openrouter, gemini-flash-lite,
 * groq (gpt-oss-20b), claude-opus-4-8, claude-sonnet-5, claude-haiku-4-5
 *
 * Usage: node skills_shim.js --payload-file <path>
 * Input JSON: { model: string, prompt: string, systemContext?: string }
 *
 * HOW THIS DIFFERS FROM baseline_shim.js / recipe_benchmark_runner.js
 * ----------------------------------------------------------------------
 * The earlier recipe-injection condition pasted the ENTIRE recipe.md text
 * into the system prompt up front, then gave a narrow `read_recipe_file`
 * tool constrained to an enum of filenames the model was already told about
 * by name. Real Skills work differently:
 *   - Real Skills show only a name + short description up front (cheap
 *     "is this relevant?" metadata) — the full main instruction file is
 *     read on demand, via a general file-reading tool, only once selected.
 *   - Real Skills support browsing (list a directory, discover files whose
 *     names weren't given in advance) rather than only fetching a
 *     pre-announced filename.
 *
 * This shim fixes both gaps:
 *   - `systemContext` here is expected to be a SHORT metadata teaser (name +
 *     what-the-skill-does), not the full recipe.md text. Building that
 *     teaser is skills_benchmark_runner.js's job, not this file's.
 *   - The one tool exposed is `view(path)` — list a directory OR read a
 *     file, sandboxed to the recipes/ folder root but otherwise
 *     unrestricted, mirroring the actual `view` tool's dual behavior. The
 *     model must discover recipe.md itself (e.g. view(".") then
 *     view("recipe.md")) before it has any procedural instructions at all,
 *     and can browse code/ the same way before reading a specific script.
 *
 * The one deliberate simplification versus real Skills: SKILL SELECTION
 * ITSELF is still skipped — this condition is still "oracle-selected", i.e.
 * the teaser text tells the model this skill IS the relevant one, rather
 * than making it choose among several candidate skills. Modeling selection
 * among multiple skills is a further, separate extension, not built here.
 *
 * IMPORTANT DISCLOSURE — the bounded-exploration-then-forced-final-answer
 * mechanism below (MAX_TOOL_TURNS, forceFinalOpenAIAnswer(),
 * forceFinalAnthropicAnswer()) is NOT documented Claude Skills behavior.
 * Anthropic's tool-use docs describe an ordinary request/respond loop with
 * no arbitrary turn cap or mandatory tools-disabled final turn. This is
 * THIS BENCHMARK'S OWN cost/robustness policy — added so a model that's
 * genuinely still making forward progress (reading new files, not looping)
 * isn't discarded as a hard failure purely because it needed more turns
 * than an arbitrary cap allows. Every attempt that used it is flagged via
 * `tool_budget_exhausted` so it stays distinguishable from a model that
 * read a few relevant files and answered voluntarily — do not conflate the
 * two when interpreting results.
 *
 * CRITICAL: exactly like baseline_shim.js's tool, `view` is ONLY attached
 * when systemContext is present. If it were always available, the vanilla
 * (no-skill) baseline would silently gain filesystem access to a skill it's
 * never told exists — that would contaminate the control condition.
 *
 * Output JSON: { result, error, alreadyExtracted, raw_preview, tool_calls_made }
 * tool_calls_made here logs the actual {action, path} pairs requested, in
 * order, e.g. [{action:"list",path:"."}, {action:"read",path:"recipe.md"}, ...]
 */

import "dotenv/config";
import fs   from "fs";
import path from "path";
import { Buffer } from "buffer";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ── helpers ──────────────────────────────────────────────────────────────────

const SYSTEM_PROMPT = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.`;

const SYSTEM_PROMPT_WITH_TOOLS_SUFFIX = `

You have access to a "view" tool for exploring a relevant skill's folder. Call view(path) with path="." (or any directory) to list its contents, or with a file's path to read its full contents. You do not yet have the skill's full instructions — use view to find and read its main instruction file first. Then inspect only the smallest number of supporting files needed for the SPECIFIC task you were given — usually one to three. Do not read every pipeline script the recipe mentions; only read what this particular task actually requires. Once you have enough information, stop calling tools. Your FINAL response must still be ONLY a single fenced R code block as described above — no extra text, and no further tool calls after you start writing the final answer.`;

function combineSystem(systemContext, toolsEnabled) {
  if (!systemContext) return SYSTEM_PROMPT;
  return `${SYSTEM_PROMPT}\n\n${systemContext}${toolsEnabled ? SYSTEM_PROMPT_WITH_TOOLS_SUFFIX : ""}`;
}

// ── sandboxed general file/directory access (the "view" tool) ───────────────

// skills_shim.js lives at benchmark/scripts/ — recipes/ is a sibling of
// benchmark/, same convention used throughout this project.
const RECIPES_ROOT = path.resolve(__dirname, "..", "..", "recipes");
// path.resolve() above only does LEXICAL normalisation (collapsing ".."
// segments, joining) — it never touches the filesystem, so it cannot detect
// a symlink inside recipes/ that points somewhere else entirely. Resolve
// the real root once at startup so every request can be checked against
// where recipes/ ACTUALLY is on disk, not just its lexical path.
let RECIPES_ROOT_REAL;
try {
  RECIPES_ROOT_REAL = fs.realpathSync(RECIPES_ROOT);
} catch {
  RECIPES_ROOT_REAL = RECIPES_ROOT; // recipes/ itself missing — viewPath()'s own stat call will report this clearly
}
const MAX_FILE_BYTES = 100_000; // guard against one call blowing the context budget
const MAX_LIST_ENTRIES = 200;

// Resolves `relPath` against RECIPES_ROOT and refuses anything that would
// escape it — this is the ONLY safety boundary now (no enum, since the
// whole point is unrestricted browsing within the folder). ".." segments
// and absolute paths are caught by the lexical containment check below.
// Symlink escapes are caught SEPARATELY, by resolving the REAL path (via
// fs.realpathSync, which does follow symlinks) and re-checking containment
// against the real root — a symlink inside recipes/ pointing elsewhere
// would pass the lexical check but get caught here.
function resolveSandboxed(relPath) {
  const target = relPath && relPath.trim() ? relPath.trim() : ".";
  const fullPath = path.resolve(RECIPES_ROOT, target);
  if (fullPath !== RECIPES_ROOT && !fullPath.startsWith(RECIPES_ROOT + path.sep)) {
    return { error: `Path escapes the skill folder: "${relPath}"` };
  }

  let realFullPath;
  try {
    realFullPath = fs.realpathSync(fullPath);
  } catch {
    // Doesn't exist (yet) — let the caller's stat/read produce a clean
    // "nothing found" error rather than a confusing symlink-check failure.
    return { fullPath };
  }
  if (realFullPath !== RECIPES_ROOT_REAL && !realFullPath.startsWith(RECIPES_ROOT_REAL + path.sep)) {
    return { error: `Path escapes the skill folder via a symlink: "${relPath}"` };
  }
  return { fullPath: realFullPath };
}

function looksBinary(buf) {
  // Cheap heuristic: a NUL byte in the first few KB almost always means
  // binary content. Good enough to avoid returning garbled output for
  // anything unexpected in the folder; this benchmark's actual files are
  // all plain text (.R, .json, .txt, .md).
  const sample = buf.subarray(0, Math.min(buf.length, 8000));
  return sample.includes(0);
}

function viewPath(relPath) {
  const resolved = resolveSandboxed(relPath);
  if (resolved.error) return resolved;
  const { fullPath } = resolved;

  let stat;
  try {
    stat = fs.statSync(fullPath);
  } catch (e) {
    return { error: `Nothing found at "${relPath}": ${e.message}` };
  }

  if (stat.isDirectory()) {
    let entries;
    try {
      entries = fs.readdirSync(fullPath, { withFileTypes: true });
    } catch (e) {
      return { error: `Could not list directory "${relPath}": ${e.message}` };
    }
    const listing = entries
      .slice(0, MAX_LIST_ENTRIES)
      .map(e => `${e.isDirectory() ? "[dir] " : "[file]"} ${e.name}`)
      .sort()
      .join("\n");
    return { content: `Directory listing for "${relPath}":\n${listing || "(empty)"}` };
  }

  // File
  let buf;
  try {
    buf = fs.readFileSync(fullPath);
  } catch (e) {
    return { error: `Could not read "${relPath}": ${e.message}` };
  }
  if (looksBinary(buf)) {
    return { error: `"${relPath}" appears to be a binary file — not readable as text.` };
  }
  let text = buf.toString("utf8");
  let truncatedNote = "";
  if (buf.length > MAX_FILE_BYTES) {
    text = buf.subarray(0, MAX_FILE_BYTES).toString("utf8");
    truncatedNote = `\n\n[TRUNCATED — file is ${buf.length} bytes, showing first ${MAX_FILE_BYTES}]`;
  }
  return { content: text + truncatedNote };
}

const VIEW_TOOL_DESCRIPTION =
  "View a file's contents or list a directory's contents within the current skill's folder. Pass a directory path (e.g. '.', 'code') to list what's inside it, or a file path (e.g. 'recipe.md', 'code/05_doublets.R') to read its full contents. Use this to discover and read the skill's main instruction file and any supporting files you need.";

const ANTHROPIC_TOOLS = [{
  name: "view",
  description: VIEW_TOOL_DESCRIPTION,
  input_schema: {
    type: "object",
    properties: {
      path: { type: "string", description: "Path relative to the skill folder root, e.g. '.', 'code', 'code/05_doublets.R', 'recipe.md'." },
    },
    required: ["path"],
  },
}];

const OPENAI_TOOLS = [{
  type: "function",
  function: {
    name: "view",
    description: VIEW_TOOL_DESCRIPTION,
    parameters: {
      type: "object",
      properties: {
        path: { type: "string", description: "Path relative to the skill folder root, e.g. '.', 'code', 'code/05_doublets.R', 'recipe.md'." },
      },
      required: ["path"],
    },
  },
}];

// was a fixed 8 — raised default to 12 and made configurable after real k=1
// data (gemini-openrouter on medium_02/hard_01) showed a model genuinely
// needing more than 8 turns of legitimate, non-repeated exploration before
// it would have been ready to answer. Still bounded, not unlimited — see
// the forced-final-synthesis mechanism below, which is the real fix for
// "ran out of budget" vs. raising this number alone (which a sufficiently
// thorough model could still exceed).
const MAX_TOOL_TURNS = Number(process.env.BENCH_MAX_TOOL_TURNS || 12);

// ── retry on transient errors (identical policy to baseline_shim.js) ────────
async function withRetry(fn, { retries = 5, baseDelayMs = 2000 } = {}) {
  let lastErr;
  for (let attempt = 0; attempt <= retries; attempt++) {
    try {
      return await fn();
    } catch (e) {
      lastErr = e;
      const status = e?.status || e?.response?.status;
      const message = e?.message || String(e);
      const isRateLimit  = status === 429 || /rate.?limit/i.test(message);
      const isOverloaded = status === 529 || status === 503 || /overloaded/i.test(message);
      const isRetryable  = isRateLimit || isOverloaded;
      if (!isRetryable || attempt === retries) throw e;

      let waitMs = baseDelayMs * Math.pow(2, attempt);
      const msMatch = message.match(/try again in ([\d.]+)\s*ms/i);
      const sMatch  = message.match(/try again in ([\d.]+)\s*s(?:ec)?\b/i);
      if (msMatch)      waitMs = Math.ceil(parseFloat(msMatch[1])) + 250;
      else if (sMatch)  waitMs = Math.ceil(parseFloat(sMatch[1]) * 1000) + 250;

      const reason = isRateLimit ? "rate limit" : "server overloaded";
      process.stderr.write(`[skills_shim] ${reason} hit (attempt ${attempt + 1}/${retries + 1}), retrying in ${waitMs}ms — ${message.slice(0, 200)}\n`);
      await new Promise(r => setTimeout(r, waitMs));
    }
  }
  throw lastErr;
}

function extractRCode(text) {
  if (!text) return null;
  text = text.replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "").trim();
  let match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  match = text.match(/```[a-z]*\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) return text;
  return null;
}

// Attaches whatever tool calls were already accumulated before a failure to
// the Error object itself. Without this, throwing mid-loop silently
// discards the local toolCallsMade array — the caller only ever sees data
// returned via a normal `return`, not what existed in the stack frame at
// the moment of an exception. That made "exceeded max tool turns" errors
// look like zero exploration happened (tool_calls_made: []) when in fact
// hitting the turn cap requires having called tools on every prior turn —
// real, extensive exploration that we'd otherwise lose all visibility into
// right when it's most useful to understand what happened.
function throwWithToolCalls(message, toolCallsMade) {
  const err = new Error(message);
  err.toolCallsMade = toolCallsMade;
  throw err;
}

// Forces a final answer once the exploration budget is exhausted, instead
// of throwing away a model that was genuinely making progress (reading new
// files each turn, never repeating) and just needed more room than the cap
// allows. Deliberately omits `tools` entirely rather than relying on
// `tool_choice: "none"` — this function serves 5 different
// OpenAI-compatible endpoints (OpenAI x2, OpenRouter x2, Groq) and not all
// of them are guaranteed to handle tool_choice semantics identically;
// simply not offering any tools at all is unambiguous everywhere.
async function forceFinalOpenAIAnswer({ model, client, messages, extraParams, toolCallsMade }) {
  messages.push({
    role: "user",
    content: "The skill-file exploration budget is now exhausted. Do not request any more tools. Using the instructions and files already retrieved, produce the final executable R code now.",
  });

  const response = await client.chat.completions.create({
    model,
    messages,
    ...extraParams,
  });

  const choice = response.choices?.[0];
  if (!choice) {
    throwWithToolCalls(`${model} returned no choice during forced final synthesis`, toolCallsMade);
  }
  if (choice.finish_reason === "length") {
    throwWithToolCalls(`${model} hit the token limit during forced final synthesis — usage: ${JSON.stringify(response.usage)}`, toolCallsMade);
  }
  return { text: choice.message?.content, toolCallsMade, tool_budget_exhausted: true };
}

// ── shared OpenAI-compatible tool-calling loop ───────────────────────────────
async function callOpenAICompatibleWithTools({ model, client, prompt, systemContext, extraParams = {} }) {
  const toolsEnabled = !!systemContext;
  const messages = [
    { role: "system", content: combineSystem(systemContext, toolsEnabled) },
    { role: "user",   content: prompt },
  ];
  const toolCallsMade = [];
  let forceFinalNextTurn = false;

  // Runs one extra iteration past MAX_TOOL_TURNS: the model gets exactly
  // MAX_TOOL_TURNS turns to explore, and — if it still wants more on the
  // last of those — one additional forced, tools-disabled turn to actually
  // answer with what it's already retrieved. This also respects the tool-
  // calling protocol correctly: every tool_call the model emits gets its
  // tool_result appended BEFORE the next completion request, even on the
  // final allowed turn — the forced call happens strictly after that.
  for (let turn = 0; turn <= MAX_TOOL_TURNS + 1; turn++) {
    if (forceFinalNextTurn) {
      return forceFinalOpenAIAnswer({ model, client, messages, extraParams, toolCallsMade });
    }

    const response = await client.chat.completions.create({
      model,
      messages,
      ...(toolsEnabled ? { tools: OPENAI_TOOLS } : {}),
      ...extraParams,
    });
    const choice = response.choices[0];
    const msg = choice.message;

    if (msg.tool_calls && msg.tool_calls.length > 0) {
      messages.push(msg);
      for (const call of msg.tool_calls) {
        let relPath = null;
        try { relPath = JSON.parse(call.function.arguments).path; } catch { /* handled by viewPath's own error */ }
        const { content, error } = viewPath(relPath);
        toolCallsMade.push({ path: relPath, result: error ? "error" : "ok" });
        messages.push({
          role: "tool",
          tool_call_id: call.id,
          content: error ? `ERROR: ${error}` : content,
        });
      }
      if (turn >= MAX_TOOL_TURNS - 1) forceFinalNextTurn = true; // MAX_TOOL_TURNS counts exploration turns exactly (0-indexed loop, so this fires after the (MAX_TOOL_TURNS)th turn)
      continue;
    }

    if (choice.finish_reason === "length") {
      throwWithToolCalls(`${model} hit token limit before completing the answer — usage: ${JSON.stringify(response.usage)}`, toolCallsMade);
    }
    return { text: msg.content, toolCallsMade, tool_budget_exhausted: false };
  }
}

// ── providers ─────────────────────────────────────────────────────────────────

async function callGPT55(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  return callOpenAICompatibleWithTools({
    model: "gpt-5.5", client, prompt, systemContext,
    extraParams: { max_completion_tokens: 16384 },
  });
}

async function callGPT54Mini(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  return callOpenAICompatibleWithTools({
    model: "gpt-5.4-mini", client, prompt, systemContext,
    extraParams: { temperature: 0.6, max_completion_tokens: 8192 },
  });
}

async function callGeminiOpenRouter(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENROUTER_API_KEY, baseURL: "https://openrouter.ai/api/v1" });
  return callOpenAICompatibleWithTools({
    model: "google/gemini-3.5-flash", client, prompt, systemContext,
    extraParams: { temperature: 0.6, max_tokens: 8192 },
  });
}

async function callGeminiFlashLite(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENROUTER_API_KEY, baseURL: "https://openrouter.ai/api/v1" });
  return callOpenAICompatibleWithTools({
    model: "google/gemini-3.1-flash-lite", client, prompt, systemContext,
    extraParams: { temperature: 0.6, max_tokens: 8192 },
  });
}

async function callGroq(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.GROQ_API_KEY, baseURL: "https://api.groq.com/openai/v1" });
  return callOpenAICompatibleWithTools({
    model: "openai/gpt-oss-20b", client, prompt, systemContext,
    // Kept modest, same rationale as baseline_shim.js: Groq's free tier is
    // 8000 TPM total, and true progressive disclosure means MORE tool turns
    // per attempt on average than the narrower recipe-injection condition
    // (list, then read metadata, then maybe list code/, then read a
    // script...), which eats into that budget turn by turn. withRetry()
    // handles resulting 429s but this is a real, documented constraint.
    extraParams: { temperature: 0.6, max_tokens: 3072 },
  });
}

// Anthropic equivalent of forceFinalOpenAIAnswer — same rationale: give a
// model that was genuinely making progress one final, tools-disabled turn
// to answer with what it already retrieved, instead of discarding it.
// combineSystem(..., false) is used deliberately here (not the toolsEnabled
// value from the outer call) — the system prompt's "you have a view tool"
// suffix would be actively misleading on a turn where tools are no longer
// offered.
async function forceFinalAnthropicAnswer({ model, client, systemContext, messages, toolCallsMade }) {
  messages.push({
    role: "user",
    content: "The skill-file exploration budget is now exhausted. Do not request any more tools. Using the instructions and files already retrieved, produce the final executable R code now.",
  });

  const response = await client.messages.create({
    model,
    max_tokens: 8192,
    system: combineSystem(systemContext, false),
    messages,
  });

  if (response.stop_reason === "max_tokens") {
    throwWithToolCalls(`${model} hit the token limit during forced final synthesis (stop_reason=max_tokens). usage: ${JSON.stringify(response.usage)}`, toolCallsMade);
  }
  const textBlocks = (response.content || []).filter(b => b.type === "text").map(b => b.text);
  const text = textBlocks.join("\n").trim();
  if (!text) {
    const blockTypes = (response.content || []).map(b => b.type).join(", ") || "(no content blocks)";
    throwWithToolCalls(`${model} returned no text content during forced final synthesis (stop_reason=${response.stop_reason}, content block types: [${blockTypes}])`, toolCallsMade);
  }
  return { text, toolCallsMade, tool_budget_exhausted: true };
}

async function callAnthropicClaude(model, prompt, systemContext) {
  const Anthropic = (await import("@anthropic-ai/sdk")).default;
  const client = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });
  const toolsEnabled = !!systemContext;

  const messages = [{ role: "user", content: prompt }];
  const toolCallsMade = [];
  let forceFinalNextTurn = false;

  for (let turn = 0; turn <= MAX_TOOL_TURNS + 1; turn++) {
    if (forceFinalNextTurn) {
      return forceFinalAnthropicAnswer({ model, client, systemContext, messages, toolCallsMade });
    }

    const response = await client.messages.create({
      model,
      max_tokens: 8192,
      system: combineSystem(systemContext, toolsEnabled),
      ...(toolsEnabled ? { tools: ANTHROPIC_TOOLS } : {}),
      messages,
    });

    if (response.stop_reason === "max_tokens") {
      throwWithToolCalls(`${model} hit max_tokens before finishing its response (stop_reason=max_tokens, turn ${turn}). usage: ${JSON.stringify(response.usage)}`, toolCallsMade);
    }

    if (response.stop_reason === "tool_use") {
      messages.push({ role: "assistant", content: response.content });
      const toolResults = [];
      for (const block of response.content) {
        if (block.type !== "tool_use") continue;
        const relPath = block.input?.path;
        const { content, error } = viewPath(relPath);
        toolCallsMade.push({ path: relPath, result: error ? "error" : "ok" });
        toolResults.push({
          type: "tool_result",
          tool_use_id: block.id,
          content: error ? `ERROR: ${error}` : content,
          is_error: !!error,
        });
      }
      messages.push({ role: "user", content: toolResults });
      if (turn >= MAX_TOOL_TURNS - 1) forceFinalNextTurn = true; // MAX_TOOL_TURNS counts exploration turns exactly (0-indexed loop, so this fires after the (MAX_TOOL_TURNS)th turn)
      continue;
    }

    const textBlocks = (response.content || []).filter(b => b.type === "text").map(b => b.text);
    const text = textBlocks.join("\n").trim();
    if (!text) {
      const blockTypes = (response.content || []).map(b => b.type).join(", ") || "(no content blocks)";
      throwWithToolCalls(`${model} returned no text content (stop_reason=${response.stop_reason}, content block types: [${blockTypes}])`, toolCallsMade);
    }
    return { text, toolCallsMade, tool_budget_exhausted: false };
  }
}

// ── main ──────────────────────────────────────────────────────────────────────

function loadPayload() {
  const args = process.argv.slice(2);
  if (args[0] === "--payload-file") {
    return JSON.parse(fs.readFileSync(args[1], "utf8"));
  }
  return JSON.parse(Buffer.from(args[0], "base64").toString("utf8"));
}

const payload = loadPayload();
const { model, prompt, systemContext } = payload;

if (systemContext) {
  process.stderr.write(`[skills_shim:${model}] systemContext (teaser) present (${systemContext.length} chars) — view tool enabled for this call\n`);
}

let result = null;
let error  = null;
let raw    = null;
let toolCallsMade = [];

let toolBudgetExhausted = false;

try {
  let callResult;
  if      (model === "gpt-5.5")            { callResult = await withRetry(() => callGPT55(prompt, systemContext)); }
  else if (model === "gpt-5.4-mini")       { callResult = await withRetry(() => callGPT54Mini(prompt, systemContext)); }
  else if (model === "gemini-openrouter")   { callResult = await withRetry(() => callGeminiOpenRouter(prompt, systemContext)); }
  else if (model === "gemini-flash-lite")   { callResult = await withRetry(() => callGeminiFlashLite(prompt, systemContext)); }
  else if (model === "claude-opus-4-8")     { callResult = await withRetry(() => callAnthropicClaude("claude-opus-4-8", prompt, systemContext)); }
  else if (model === "claude-sonnet-5")     { callResult = await withRetry(() => callAnthropicClaude("claude-sonnet-5", prompt, systemContext)); }
  else if (model === "claude-haiku-4-5")    { callResult = await withRetry(() => callAnthropicClaude("claude-haiku-4-5-20251001", prompt, systemContext)); }
  else if (model === "groq")                { callResult = await withRetry(() => callGroq(prompt, systemContext)); }
  else throw new Error(`Unknown or unsupported model for skills_shim.js: ${model}`);

  raw = callResult.text;
  toolCallsMade = callResult.toolCallsMade || [];
  toolBudgetExhausted = !!callResult.tool_budget_exhausted;
  if (toolCallsMade.length > 0) {
    process.stderr.write(`[skills_shim:${model}] view() calls: ${toolCallsMade.map(t => `${t.path}(${t.result})`).join(", ")}${toolBudgetExhausted ? " [forced final synthesis after budget exhausted]" : ""}\n`);
  }

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in response`);

} catch (e) {
  error = e.message;
  // If the error was thrown mid-loop (see throwWithToolCalls), it carries
  // whatever exploration already happened before the failure — surface
  // that instead of leaving toolCallsMade at its initial empty array,
  // otherwise a genuine multi-turn exploration that ultimately failed
  // looks identical to "never explored at all" in the saved results.
  if (Array.isArray(e.toolCallsMade)) {
    toolCallsMade = e.toolCallsMade;
  }
  process.stderr.write(`[skills_shim:${model}] Error: ${error}\n`);
  if (toolCallsMade.length > 0) {
    process.stderr.write(`[skills_shim:${model}] view() calls before failure: ${toolCallsMade.map(t => `${t.path}(${t.result})`).join(", ")}\n`);
  }
  if (raw) {
    process.stderr.write(`[skills_shim:${model}] Raw response (first 1000 chars):\n${String(raw).slice(0, 1000)}\n`);
  }
}

process.stdout.write(JSON.stringify({
  result,
  error,
  alreadyExtracted: true,
  raw_preview: (!result && raw) ? String(raw).slice(0, 500) : undefined,
  tool_calls_made: toolCallsMade,
  tool_budget_exhausted: toolBudgetExhausted,
}) + "\n");
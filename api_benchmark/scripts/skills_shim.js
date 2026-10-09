// skills_shim.js
// Like baseline_shim.js, plus the view() tool loop for the "With Skills" condition.
// Prints one JSON line: { result, error, alreadyExtracted, raw_preview,
// tool_calls_made, tool_budget_exhausted }.
//
// Usage: node skills_shim.js --payload-file <path>
// Input: { model, prompt, systemContext? }   (systemContext = the short skill teaser)
//
// view(path) lists a directory or reads a file, restricted to recipes/. The tool is
// only offered when systemContext is present, so the vanilla condition never has it.
// The model gets up to MAX_TOOL_TURNS turns to explore and then one final turn
// without tools; such attempts are flagged with tool_budget_exhausted.

import "dotenv/config";
import fs   from "fs";
import path from "path";
import { Buffer } from "buffer";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ---- helpers ----

const SYSTEM_PROMPT = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.`;

const SYSTEM_PROMPT_WITH_TOOLS_SUFFIX = `

You have access to a "view" tool for exploring a relevant skill's folder. Call view(path) with path="." (or any directory) to list its contents, or with a file's path to read its full contents. You do not yet have the skill's full instructions — use view to find and read its main instruction file first. Then inspect only the smallest number of supporting files needed for the SPECIFIC task you were given — usually one to three. Do not read every pipeline script the recipe mentions; only read what this particular task actually requires. Once you have enough information, stop calling tools. Your FINAL response must still be ONLY a single fenced R code block as described above — no extra text, and no further tool calls after you start writing the final answer.`;

// system prompt + teaser, plus the view-tool note when tools are offered on this turn
function combineSystem(systemContext, toolsEnabled) {
  if (!systemContext) return SYSTEM_PROMPT;
  return `${SYSTEM_PROMPT}\n\n${systemContext}${toolsEnabled ? SYSTEM_PROMPT_WITH_TOOLS_SUFFIX : ""}`;
}

// ---- view tool (sandboxed to recipes/) ----

// recipes/ sits next to api_benchmark/
const RECIPES_ROOT = path.resolve(__dirname, "..", "..", "recipes");
// resolve the real path once so symlinks out of recipes/ can be detected
let RECIPES_ROOT_REAL;
try {
  RECIPES_ROOT_REAL = fs.realpathSync(RECIPES_ROOT);
} catch {
  RECIPES_ROOT_REAL = RECIPES_ROOT; // recipes/ missing: viewPath() will report it
}
const MAX_FILE_BYTES = 100_000; // guard against one call blowing the context budget
const MAX_LIST_ENTRIES = 200;

// Resolves relPath inside recipes/ and rejects anything that escapes it, by "..",
// absolute path or symlink.
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
    // does not exist: let the caller report "nothing found"
    return { fullPath };
  }
  if (realFullPath !== RECIPES_ROOT_REAL && !realFullPath.startsWith(RECIPES_ROOT_REAL + path.sep)) {
    return { error: `Path escapes the skill folder via a symlink: "${relPath}"` };
  }
  return { fullPath: realFullPath };
}

function looksBinary(buf) {
  // a NUL byte in the first few KB means binary
  const sample = buf.subarray(0, Math.min(buf.length, 8000));
  return sample.includes(0);
}

// The view tool: directory listing or file contents, unmodified (truncated at MAX_FILE_BYTES).
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

// raised from 8: some models needed more turns of real exploration
const MAX_TOOL_TURNS = Number(process.env.BENCH_MAX_TOOL_TURNS || 12);

// ---- retry on transient errors (same as baseline_shim.js) ----
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

// Throws an error that carries the view() calls made so far, so a failed attempt
// still records its exploration.
function throwWithToolCalls(message, toolCallsMade) {
  const err = new Error(message);
  err.toolCallsMade = toolCallsMade;
  throw err;
}

// Tool budget used up: ask for a final answer with no tools offered at all
// (not every endpoint handles tool_choice: "none" the same way).
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

// ---- OpenAI-compatible tool loop ----
// Plain completion without systemContext. With it, the model may call view() for up to
// MAX_TOOL_TURNS turns, then gets one forced final turn.
async function callOpenAICompatibleWithTools({ model, client, prompt, systemContext, extraParams = {} }) {
  const toolsEnabled = !!systemContext;
  const messages = [
    { role: "system", content: combineSystem(systemContext, toolsEnabled) },
    { role: "user",   content: prompt },
  ];
  const toolCallsMade = [];
  let forceFinalNextTurn = false;

  // one extra iteration: every tool call gets its result before the forced final turn
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
      if (turn >= MAX_TOOL_TURNS - 1) forceFinalNextTurn = true; // last exploration turn
      continue;
    }

    if (choice.finish_reason === "length") {
      throwWithToolCalls(`${model} hit token limit before completing the answer — usage: ${JSON.stringify(response.usage)}`, toolCallsMade);
    }
    return { text: msg.content, toolCallsMade, tool_budget_exhausted: false };
  }
}

// ---- providers ----

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

// Same as forceFinalOpenAIAnswer for Anthropic. The system prompt is built without the
// view-tool note because no tools are offered on this turn.
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

// Anthropic tool loop, same policy as the OpenAI one.
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
      if (turn >= MAX_TOOL_TURNS - 1) forceFinalNextTurn = true; // last exploration turn
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

// ---- main ----

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
  // errors thrown mid-loop carry the view() calls made before the failure
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
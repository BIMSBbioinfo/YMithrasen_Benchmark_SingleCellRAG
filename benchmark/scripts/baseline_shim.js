/**
 * baseline_shim.js
 *
 * Calls various model providers and returns extracted R code.
 * Supports: gpt-5.5, gpt-5.4-mini, gemini-openrouter, gemini-flash-lite,
 * groq (gpt-oss-20b), claude-opus-4-8, claude-sonnet-5, claude-haiku-4-5,
 * rolv-raw
 *
 * Usage: node baseline_shim.js <base64-encoded-json>
 * Input JSON: { model: string, prompt: string, systemContext?: string }
 *
 * systemContext (optional): additional context — e.g. a recipe/skill
 * document — injected into the provider's system/developer instruction
 * channel, NOT concatenated into the user prompt. When present, it is
 * appended after the fixed SYSTEM_PROMPT (below) rather than replacing it,
 * so the output-format contract ("respond with ONLY executable R code,
 * fenced") stays identical between the vanilla and recipe-injected
 * conditions — the only thing that changes is the presence of the recipe.
 * Routing per provider:
 *   - Anthropic (claude-*):        `system` request param
 *   - OpenAI-compatible (GPT,
 *     Gemini via OpenRouter, Groq): a `role: "system"` message
 *   - rolv-raw:                    a `role: "system"` message (same as the
 *                                  OpenAI-compatible providers, since the
 *                                  ROLV raw endpoint is OpenAI-chat-shaped)
 *   - rolv (grounded/RAG mode):    intentionally NOT wired to systemContext
 *                                  — see callRolv() below; out of scope for
 *                                  the recipe-injection benchmark anyway.
 *
 * Sampling note: temperature is set to 0.6 for providers that support it
 * (Gemini, Groq, GPT-5.4-mini) to give genuine k-sample diversity for Pass@k.
 * GPT-5.5 and Claude direct models do not accept temperature — they run at
 * provider default. This asymmetry is a real, documented limitation: not all
 * models received identical sampling controls. Diversity across k attempts
 * is tracked separately (see unique_generations in benchmark_runner.js).
 * Output JSON: { result: string|null, error: string|null }
 */

import "dotenv/config";
import fs from "fs";
import { Buffer } from "buffer";

// ── helpers ──────────────────────────────────────────────────────────────────

const SYSTEM_PROMPT = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.`;

// Appends systemContext (e.g. a recipe document) after the fixed
// SYSTEM_PROMPT rather than replacing it, so the output-format contract is
// identical between vanilla and recipe-injected runs — only the presence of
// the recipe changes. Returns SYSTEM_PROMPT unchanged when systemContext is
// absent, so vanilla-baseline calls are byte-identical to before this patch.
function combineSystem(systemContext) {
  return systemContext ? `${SYSTEM_PROMPT}\n\n${systemContext}` : SYSTEM_PROMPT;
}

// Retries a provider call on transient errors instead of letting them fall
// straight through as a scored failure:
//   - 429 rate limits — systemContext added ~2-3K tokens to EVERY call,
//     which can push providers with tight per-minute budgets (Groq's free
//     tier: 8000 TPM) over their limit within a couple of calls even though
//     the model itself never got a real chance to answer.
//   - 5xx / "overloaded" server errors (e.g. Anthropic's 529
//     "overloaded_error") — transient capacity issues on the provider's
//     side, unrelated to the prompt, the recipe, or the model's actual
//     capability. A single un-retried batch of these can silently wipe out
//     a large fraction of one model's prompts and make it look like a real
//     capability gap when it's actually "the API was down for 40 seconds".
// Scoring either as a Pass@k failure would corrupt the benchmark for
// reasons that have nothing to do with recipe injection or the model.
// Parses the provider's suggested wait time out of the error message when
// present (e.g. Groq's "Please try again in 794.999999ms"); falls back to
// exponential backoff otherwise (used for 529s, which typically don't state
// a wait time). Non-retryable errors are rethrown immediately.
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
      if (msMatch)      waitMs = Math.ceil(parseFloat(msMatch[1])) + 250; // small buffer
      else if (sMatch)  waitMs = Math.ceil(parseFloat(sMatch[1]) * 1000) + 250;

      const reason = isRateLimit ? "rate limit" : "server overloaded";
      process.stderr.write(`[baseline_shim] ${reason} hit (attempt ${attempt + 1}/${retries + 1}), retrying in ${waitMs}ms — ${message.slice(0, 200)}\n`);
      await new Promise(r => setTimeout(r, waitMs));
    }
  }
  throw lastErr;
}

function extractRCode(text) {
  if (!text) return null;
  // Strip scratchpad blocks
  text = text.replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "").trim();
  // Try R-specific fence first
  let match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  // Fall back to any code fence
  match = text.match(/```[a-z]*\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  // If no fence, treat whole text as R code if it looks like R
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) {
    return text;
  }
  return null;
}

// ── providers ─────────────────────────────────────────────────────────────────

async function callGPT55(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-5.5",
    max_completion_tokens: 16384, // was 2048 — reasoning tokens were consuming the entire budget on complex prompts, leaving 0 tokens for the visible answer
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  const choice = response.choices[0];
  if (!choice.message.content && choice.finish_reason === "length") {
    throw new Error(`GPT-5.5 hit max_completion_tokens before producing visible output (reasoning tokens likely consumed the budget) — usage: ${JSON.stringify(response.usage)}`);
  }
  return choice.message.content;
}

async function callGPT54Mini(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-5.4-mini",
    temperature: 0.6,
    max_completion_tokens: 8192, // was 2048 — same reasoning-token headroom issue as GPT-5.5
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGeminiOpenRouter(prompt, systemContext) {
  // Gemini 3.5 Flash via OpenRouter
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.OPENROUTER_API_KEY,
    baseURL: "https://openrouter.ai/api/v1"
  });
  const response = await client.chat.completions.create({
    model: "google/gemini-3.5-flash",
    temperature: 0.6,
    max_tokens: 8192, // was 2048 — same token-starvation bug already fixed for GPT-5.5/5.4-mini: recipe-injected prompts push the model toward longer, more complete pipelines (data load -> QC -> doublets -> normalization -> PCA -> UMAP -> clustering all in one response), which routinely blew past 2048 and got cut off mid-line. Confirmed via recipe_results_2026-07-29T17-24-39-162.json: 8/30 gemini-openrouter responses truncated before a closing code fence.
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGeminiFlashLite(prompt, systemContext) {
  // Gemini 3.1 Flash-Lite via OpenRouter — cheapest Gemini, GA version
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.OPENROUTER_API_KEY,
    baseURL: "https://openrouter.ai/api/v1"
  });
  const response = await client.chat.completions.create({
    model: "google/gemini-3.1-flash-lite",
    temperature: 0.6,
    max_tokens: 8192, // was 2048 — see callGeminiOpenRouter comment; raised preemptively for the same reason even though Flash-Lite showed 0 truncations in the k=1 dry run (it tends to write terser code, but there's no guarantee that holds across all 29 prompts at k=3)
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callAnthropicClaude(model, prompt, systemContext) {
  // Direct Anthropic API — bypasses ROLV entirely, uses ANTHROPIC_API_KEY
  // Note: Sonnet 5, Opus 4.7, Opus 4.8 do NOT support temperature/top_p/top_k
  const Anthropic = (await import("@anthropic-ai/sdk")).default;
  const client = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });
  const response = await client.messages.create({
    model,
    max_tokens: 8192, // was 2048 — same token-starvation bug already fixed for GPT-5.5/5.4-mini, never applied here. Confirmed via recipe_results_2026-07-29T17-24-39-162.json: 13/90 Claude-direct responses (across all 3 models) were cut off mid-line before a closing code fence, because recipe injection pushes toward longer, more complete pipelines that routinely exceed 2048 tokens.
    system: combineSystem(systemContext),
    messages: [{ role: "user", content: prompt }]
  });
  if (response.stop_reason === "max_tokens") {
    throw new Error(`${model} hit max_tokens before finishing its response (stop_reason=max_tokens) — output was likely truncated mid-code. usage: ${JSON.stringify(response.usage)}`);
  }
  return response.content[0].text;
}

async function callGroq(prompt, systemContext) {
  // Groq — OpenAI GPT-OSS 20B (free at console.groq.com/keys)
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.GROQ_API_KEY,
    baseURL: "https://api.groq.com/openai/v1"
  });
  const response = await client.chat.completions.create({
    model: "openai/gpt-oss-20b",
    temperature: 0.6,
    max_tokens: 3072, // was 2048. NOT raised as far as the other providers (8192) deliberately: Groq's free tier caps at 8000 TPM total (prompt + completion), and systemContext alone already consumes ~2700 of that per call — requesting 8192 completion tokens would frequently exceed the entire per-minute budget in a single call. 3072 is a compromise: reduces truncation risk somewhat without making every call rate-limit-bound. withRetry() will retry any 429s this still causes, but a persistently tight budget is a real, documented limitation of this model under the recipe-injection condition — worth noting in the write-up rather than something further token-limit tuning can fully solve.
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callRolv(prompt) {
  // ROLV grounded mode — use ROLV's own system prompt (no override)
  // Adding our own system prompt interferes with ROLV's RAG and pre-exec planner
  // Intentionally does NOT accept/apply systemContext: this mode already does
  // its own retrieval, so an external recipe injection is out of scope here
  // (see recipe_benchmark_runner.js — rolv/rolv-nopreexec are excluded from
  // that benchmark for exactly this reason). If systemContext is passed in
  // for model === "rolv" it is silently ignored, by design.
  const apiUrl = process.env.ROLV_API_URL;
  const apiKey = process.env.ROLV_JWT_TOKEN;

  const response = await fetch(apiUrl, {
    method:  "POST",
    headers: {
      "Content-Type":  "application/json",
      "Authorization": `Bearer ${apiKey}`
    },
    body: JSON.stringify({
      model: "rolv",
      max_tokens: 2048,
      temperature: 0.6,
      messages: [
        { role: "user", content: prompt }
      ]
    })
  });

  if (!response.ok) {
    const err = await response.text();
    throw new Error(`ROLV API error ${response.status}: ${err}`);
  }

  const data = await response.json();
  return data.choices?.[0]?.message?.content
      || data.content?.[0]?.text
      || JSON.stringify(data);
}

async function callRolvRaw(prompt, systemContext) {
  // ROLV raw mode — @out-of-domain bypasses RAG, use minimal system prompt.
  // The ROLV raw endpoint is OpenAI-chat-shaped (accepts a `messages` array
  // with a system-role entry), so systemContext is routed the same way as
  // the other OpenAI-compatible providers above.
  const apiUrl = process.env.ROLV_API_URL;
  const apiKey = process.env.ROLV_JWT_TOKEN;

  const response = await fetch(apiUrl, {
    method:  "POST",
    headers: {
      "Content-Type":  "application/json",
      "Authorization": `Bearer ${apiKey}`
    },
    body: JSON.stringify({
      model: "rolv",
      max_tokens: 4096,
      temperature: 0.6,
      messages: [
        { role: "system", content: combineSystem(systemContext) },
        { role: "user",   content: `@out-of-domain ${prompt}` }
      ]
    })
  });

  if (!response.ok) {
    const err = await response.text();
    throw new Error(`ROLV raw API error ${response.status}: ${err}`);
  }

  const data = await response.json();
  return data.choices?.[0]?.message?.content
      || data.content?.[0]?.text
      || JSON.stringify(data);
}

// ── main ──────────────────────────────────────────────────────────────────────

// Two invocation styles are supported:
//   node baseline_shim.js <base64-encoded-json>          (legacy — used by
//     benchmark_runner.js, where vanilla prompts are short enough to stay
//     under Windows' ~8191-char command-line length limit)
//   node baseline_shim.js --payload-file <path>          (used by
//     recipe_benchmark_runner.js, where the recipe pushes systemContext well
//     past that limit — a base64 payload embedded in the command string gets
//     silently truncated/rejected by CreateProcess on Windows, which causes
//     exec() to return near-empty output that parses as "{}" without ever
//     throwing, masking the real failure as "no R code extracted")
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
  process.stderr.write(`[baseline_shim:${model}] systemContext present (${systemContext.length} chars) — will be appended to SYSTEM_PROMPT for this call\n`);
}

let result = null;
let error  = null;
let raw    = null; // hoisted out of try{} so it survives into the catch/output below

try {
  if      (model === "rolv")               raw = await callRolv(prompt); // systemContext intentionally not accepted here — see callRolv()
  else if (model === "rolv-raw")           raw = await withRetry(() => callRolvRaw(prompt, systemContext));
  else if (model === "gpt-5.5")            raw = await withRetry(() => callGPT55(prompt, systemContext));
  else if (model === "gpt-5.4-mini")       raw = await withRetry(() => callGPT54Mini(prompt, systemContext));
  else if (model === "gemini-openrouter")   raw = await withRetry(() => callGeminiOpenRouter(prompt, systemContext));
  else if (model === "gemini-flash-lite")   raw = await withRetry(() => callGeminiFlashLite(prompt, systemContext));
  else if (model === "claude-opus-4-8")     raw = await withRetry(() => callAnthropicClaude("claude-opus-4-8", prompt, systemContext));
  else if (model === "claude-sonnet-5")     raw = await withRetry(() => callAnthropicClaude("claude-sonnet-5", prompt, systemContext));
  else if (model === "claude-haiku-4-5")    raw = await withRetry(() => callAnthropicClaude("claude-haiku-4-5-20251001", prompt, systemContext));
  else if (model === "groq")                raw = await withRetry(() => callGroq(prompt, systemContext));
  else throw new Error(`Unknown model: ${model}`);

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in response`);

} catch (e) {
  error = e.message;
  process.stderr.write(`[baseline_shim:${model}] Error: ${error}\n`);
  if (raw) {
    // Log the actual raw text so extraction failures are diagnosable instead
    // of just "no code found" with no visibility into what the model said.
    process.stderr.write(`[baseline_shim:${model}] Raw response (first 1000 chars):\n${String(raw).slice(0, 1000)}\n`);
  }
}

// alreadyExtracted: true tells benchmark_runner.js this is clean R code with
// fences already stripped — do not re-run extractRFromText's heuristics on it.
// raw_preview is only populated on failure, so successful runs stay lean.
process.stdout.write(JSON.stringify({
  result,
  error,
  alreadyExtracted: true,
  raw_preview: (!result && raw) ? String(raw).slice(0, 500) : undefined,
}) + "\n");
// baseline_shim.js
// Calls one model once and prints the extracted R code as a JSON line:
// { result, error, alreadyExtracted, raw_preview }. Everything else goes to stderr.
//
// Usage: node baseline_shim.js <base64 JSON>   or   --payload-file <path>
// Input: { model, prompt, systemContext? }
//
// Models: gpt-5.5, gpt-5.4-mini, gemini-openrouter, gemini-flash-lite,
// claude-opus-4-8, claude-sonnet-5, claude-haiku-4-5, rolv, rolv-raw
//
// Sampling: temperature 0.6 where the provider accepts it (Gemini, GPT-5.4 Mini).
// GPT-5.5 and the Claude models take no temperature and run at the provider default.

import "dotenv/config";
import fs from "fs";
import { Buffer } from "buffer";

// ---- helpers ----

const SYSTEM_PROMPT = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.`;

// optional extra context is appended after the fixed system prompt, never replaces it
function combineSystem(systemContext) {
  return systemContext ? `${SYSTEM_PROMPT}\n\n${systemContext}` : SYSTEM_PROMPT;
}

// Retries on rate limits (429) and overload errors (503/529) so that provider hiccups
// are not scored as model failures. Uses the wait time from the error message if
// there is one, otherwise exponential backoff.
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

// R fence first, then any fence, then the whole reply if it looks like R
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
  // no fence: take the whole text if it looks like R
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) {
    return text;
  }
  return null;
}

// ---- providers ----

async function callGPT55(prompt, systemContext) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-5.5",
    max_completion_tokens: 16384, // reasoning tokens used up a 2048 budget
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
    max_completion_tokens: 8192, // same as GPT-5.5
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
    max_tokens: 8192, // 2048 truncated longer pipelines
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGeminiFlashLite(prompt, systemContext) {
  // Gemini 3.1 Flash-Lite via OpenRouter
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.OPENROUTER_API_KEY,
    baseURL: "https://openrouter.ai/api/v1"
  });
  const response = await client.chat.completions.create({
    model: "google/gemini-3.1-flash-lite",
    temperature: 0.6,
    max_tokens: 8192, // same as Gemini 3.5 Flash
    messages: [
      { role: "system", content: combineSystem(systemContext) },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callAnthropicClaude(model, prompt, systemContext) {
  // direct Anthropic API; these models do not accept temperature/top_p/top_k
  const Anthropic = (await import("@anthropic-ai/sdk")).default;
  const client = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });
  const response = await client.messages.create({
    model,
    max_tokens: 8192, // 2048 cut off responses mid-code
    system: combineSystem(systemContext),
    messages: [{ role: "user", content: prompt }]
  });
  if (response.stop_reason === "max_tokens") {
    throw new Error(`${model} hit max_tokens before finishing its response (stop_reason=max_tokens) — output was likely truncated mid-code. usage: ${JSON.stringify(response.usage)}`);
  }
  return response.content[0].text;
}

// ROLV over HTTP. Not used for the Agentic RAG runs, which go through the orchestrator shims.
async function callRolv(prompt) {
  // ROLV keeps its own system prompt here; systemContext is ignored on purpose
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

// ROLV with retrieval bypassed (model key "rolv-raw")
async function callRolvRaw(prompt, systemContext) {
  // @out-of-domain skips retrieval; the endpoint takes OpenAI-style messages
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

// ---- main ----

// Payload comes as a base64 argument, or as --payload-file <path> when it is too long
// for the Windows command line.
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
let raw    = null; // kept outside try so the catch block can log it

try {
  if      (model === "rolv")               raw = await callRolv(prompt); // no systemContext for rolv
  else if (model === "rolv-raw")           raw = await withRetry(() => callRolvRaw(prompt, systemContext));
  else if (model === "gpt-5.5")            raw = await withRetry(() => callGPT55(prompt, systemContext));
  else if (model === "gpt-5.4-mini")       raw = await withRetry(() => callGPT54Mini(prompt, systemContext));
  else if (model === "gemini-openrouter")   raw = await withRetry(() => callGeminiOpenRouter(prompt, systemContext));
  else if (model === "gemini-flash-lite")   raw = await withRetry(() => callGeminiFlashLite(prompt, systemContext));
  else if (model === "claude-opus-4-8")     raw = await withRetry(() => callAnthropicClaude("claude-opus-4-8", prompt, systemContext));
  else if (model === "claude-sonnet-5")     raw = await withRetry(() => callAnthropicClaude("claude-sonnet-5", prompt, systemContext));
  else if (model === "claude-haiku-4-5")    raw = await withRetry(() => callAnthropicClaude("claude-haiku-4-5-20251001", prompt, systemContext));
  else throw new Error(`Unknown model: ${model}`);

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in response`);

} catch (e) {
  error = e.message;
  process.stderr.write(`[baseline_shim:${model}] Error: ${error}\n`);
  if (raw) {
    // log the raw reply so extraction failures can be diagnosed
    process.stderr.write(`[baseline_shim:${model}] Raw response (first 1000 chars):\n${String(raw).slice(0, 1000)}\n`);
  }
}

// alreadyExtracted tells the runner not to extract again; raw_preview is only set on failure
process.stdout.write(JSON.stringify({
  result,
  error,
  alreadyExtracted: true,
  raw_preview: (!result && raw) ? String(raw).slice(0, 500) : undefined,
}) + "\n");
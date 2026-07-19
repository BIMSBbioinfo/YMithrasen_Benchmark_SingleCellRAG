/**
 * baseline_shim.js
 *
 * Calls various model providers and returns extracted R code.
 * Supports: gpt-5.5, gemini-openrouter (via OpenRouter), groq (gpt-oss-20b), rolv-raw
 *
 * Usage: node baseline_shim.js <base64-encoded-json>
 * Input JSON: { model: string, prompt: string }
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
import { Buffer } from "buffer";

// ── helpers ──────────────────────────────────────────────────────────────────

const SYSTEM_PROMPT = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.`;

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

async function callGPT55(prompt) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-5.5",
    max_completion_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGPT54Mini(prompt) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-5.4-mini",
    temperature: 0.6,
    max_completion_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGeminiOpenRouter(prompt) {
  // Gemini 3.5 Flash via OpenRouter
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.OPENROUTER_API_KEY,
    baseURL: "https://openrouter.ai/api/v1"
  });
  const response = await client.chat.completions.create({
    model: "google/gemini-3.5-flash",
    temperature: 0.6,
    max_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGeminiFlashLite(prompt) {
  // Gemini 3.1 Flash-Lite via OpenRouter — cheapest Gemini, GA version
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.OPENROUTER_API_KEY,
    baseURL: "https://openrouter.ai/api/v1"
  });
  const response = await client.chat.completions.create({
    model: "google/gemini-3.1-flash-lite",
    temperature: 0.6,
    max_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callAnthropicClaude(model, prompt) {
  // Direct Anthropic API — bypasses ROLV entirely, uses ANTHROPIC_API_KEY
  // Note: Sonnet 5, Opus 4.7, Opus 4.8 do NOT support temperature/top_p/top_k
  const Anthropic = (await import("@anthropic-ai/sdk")).default;
  const client = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });
  const response = await client.messages.create({
    model,
    max_tokens: 2048,
    system: SYSTEM_PROMPT,
    messages: [{ role: "user", content: prompt }]
  });
  return response.content[0].text;
}

async function callGroq(prompt) {
  // Groq — OpenAI GPT-OSS 20B (free at console.groq.com/keys)
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({
    apiKey:  process.env.GROQ_API_KEY,
    baseURL: "https://api.groq.com/openai/v1"
  });
  const response = await client.chat.completions.create({
    model: "openai/gpt-oss-20b",
    temperature: 0.6,
    max_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callRolv(prompt) {
  // ROLV grounded mode — use ROLV's own system prompt (no override)
  // Adding our own system prompt interferes with ROLV's RAG and pre-exec planner
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

async function callRolvRaw(prompt) {
  // ROLV raw mode — @out-of-domain bypasses RAG, use minimal system prompt
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
        { role: "system", content: SYSTEM_PROMPT },
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

const payload = JSON.parse(Buffer.from(process.argv[2], "base64").toString("utf8"));
const { model, prompt } = payload;

let result = null;
let error  = null;

try {
  let raw;
  if      (model === "rolv")               raw = await callRolv(prompt);
  else if (model === "rolv-raw")           raw = await callRolvRaw(prompt);
  else if (model === "gpt-5.5")            raw = await callGPT55(prompt);
  else if (model === "gpt-5.4-mini")       raw = await callGPT54Mini(prompt);
  else if (model === "gemini-openrouter")   raw = await callGeminiOpenRouter(prompt);
  else if (model === "gemini-flash-lite")   raw = await callGeminiFlashLite(prompt);
  else if (model === "claude-opus-4-8")     raw = await callAnthropicClaude("claude-opus-4-8", prompt);
  else if (model === "claude-sonnet-5")     raw = await callAnthropicClaude("claude-sonnet-5", prompt);
  else if (model === "claude-haiku-4-5")    raw = await callAnthropicClaude("claude-haiku-4-5-20251001", prompt);
  else if (model === "groq")                raw = await callGroq(prompt);
  else throw new Error(`Unknown model: ${model}`);

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in response`);

} catch (e) {
  error = e.message;
  process.stderr.write(`[baseline_shim:${model}] Error: ${error}\n`);
}

// alreadyExtracted: true tells benchmark_runner.js this is clean R code with
// fences already stripped — do not re-run extractRFromText's heuristics on it.
process.stdout.write(JSON.stringify({ result, error, alreadyExtracted: true }) + "\n");
/**
 * baseline_shim.js
 *
 * Calls various model providers and returns extracted R code.
 * Supports: gpt-4o, gemini-2.0-flash, groq (gpt-oss-20b), rolv-raw
 *
 * Usage: node baseline_shim.js <base64-encoded-json>
 * Input JSON: { model: string, prompt: string }
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

async function callGPT4o(prompt) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-4o",
    temperature: 0,
    max_tokens: 2048,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt }
    ]
  });
  return response.choices[0].message.content;
}

async function callGemini(prompt) {
  const { GoogleGenerativeAI } = await import("@google/generative-ai");
  const genAI = new GoogleGenerativeAI(process.env.GOOGLE_API_KEY);
  const model = genAI.getGenerativeModel({
    model: "gemini-2.0-flash",
    systemInstruction: SYSTEM_PROMPT,
    generationConfig: { temperature: 0, maxOutputTokens: 2048 }
  });
  const result = await model.generateContent(prompt);
  return result.response.text();
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
    temperature: 0,
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
      temperature: 0,
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
      temperature: 0,
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
  if      (model === "rolv")             raw = await callRolv(prompt);
  else if (model === "rolv-raw")         raw = await callRolvRaw(prompt);
  else if (model === "gpt-4o")           raw = await callGPT4o(prompt);
  else if (model === "gemini-2.0-flash") raw = await callGemini(prompt);
  else if (model === "groq")             raw = await callGroq(prompt);
  else throw new Error(`Unknown model: ${model}`);

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in response`);

} catch (e) {
  error = e.message;
  process.stderr.write(`[baseline_shim:${model}] Error: ${error}\n`);
}

process.stdout.write(JSON.stringify({ result, error }) + "\n");
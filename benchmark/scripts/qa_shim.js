/**
 * qa_shim.js
 *
 * Lightweight shim for Genomics QA evaluation.
 * Returns raw text responses — no R code extraction.
 * Uses Node.js built-in https — zero extra dependencies.
 *
 * payload: { model: "claude"|"groq", prompt: "..." }
 * stdout:  { result: "...", error: null }
 */

import "dotenv/config";
import https from "https";

const real = process.stdout.write.bind(process.stdout);
const err  = (msg) => process.stderr.write(msg + "\n");

// ── parse payload ─────────────────────────────────────────────────────────────

const raw = process.argv[2];
if (!raw) { err("[qa_shim] no payload"); process.exit(1); }

let model, prompt;
try {
  ({ model, prompt } = JSON.parse(Buffer.from(raw, "base64").toString("utf-8")));
} catch (e) {
  err(`[qa_shim] payload parse error: ${e.message}`);
  process.exit(1);
}

// ── http helper ───────────────────────────────────────────────────────────────

function post(hostname, path, headers, body) {
  return new Promise((resolve, reject) => {
    const data = JSON.stringify(body);
    const req  = https.request(
      {
        hostname,
        path,
        method: "POST",
        headers: { ...headers, "Content-Length": Buffer.byteLength(data) }
      },
      (res) => {
        let chunks = "";
        res.on("data", c => chunks += c);
        res.on("end", () => {
          if (res.statusCode >= 400)
            return reject(new Error(`HTTP ${res.statusCode}: ${chunks.slice(0, 200)}`));
          try { resolve(JSON.parse(chunks)); } catch { resolve(chunks); }
        });
      }
    );
    req.on("error", reject);
    req.write(data);
    req.end();
  });
}

// ── callers ───────────────────────────────────────────────────────────────────

async function callClaude(prompt) {
  // Use the ROLV proxy with model alias "rolv" — same as what baseAgent.js uses.
  // The proxy maps "rolv" → claude-sonnet-4-5 internally.
  const apiUrl = process.env.ROLV_API_URL || "https://api.rolv.io/api/chat/completions";
  const token  = process.env.ROLV_JWT_TOKEN || "";
  const url    = new URL(apiUrl);

  const data = await post(
    url.hostname,
    url.pathname,
    {
      "Content-Type":  "application/json",
      "Authorization": `Bearer ${token}`,
    },
    {
      model:      "rolv",
      messages:   [{ role: "user", content: prompt }],
      max_tokens: 512,
    }
  );
  return data?.choices?.[0]?.message?.content || "";
}

async function callGroq(prompt) {
  const apiKey = process.env.GROQ_API_KEY || "";
  const data   = await post(
    "api.groq.com",
    "/openai/v1/chat/completions",
    {
      "Content-Type":  "application/json",
      "Authorization": `Bearer ${apiKey}`,
    },
    {
      model:      "openai/gpt-oss-20b",
      messages:   [{ role: "user", content: prompt }],
      max_tokens: 512,
    }
  );
  return data?.choices?.[0]?.message?.content || "";
}

// ── main ──────────────────────────────────────────────────────────────────────

try {
  const result = (model === "groq") ? await callGroq(prompt) : await callClaude(prompt);
  real(JSON.stringify({ result, error: null }) + "\n");
} catch (e) {
  err(`[qa_shim] error: ${e.message}`);
  real(JSON.stringify({ result: null, error: e.message }) + "\n");
  process.exit(1);
}

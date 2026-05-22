/**
 * Calls baseline models (GPT-4o, Gemini 1.5 Pro) for benchmark comparison.
 */

import "dotenv/config";
import path from "path";
import { fileURLToPath } from "url";
const __dirname = path.dirname(fileURLToPath(import.meta.url));

const raw    = process.argv[2];
const parsed = JSON.parse(Buffer.from(raw, "base64").toString("utf-8"));
const { model, prompt } = parsed;

// Shared system prompt — mirror what you send to ROLV so comparison is fair
const SYSTEM_PROMPT = `You are a helpful bioinformatics assistant specializing in single-cell RNA-seq analysis using R and the Seurat package.
When asked to write R code, respond ONLY with an executable R Markdown code block using the format:
\`\`\`{r}
# your code here
\`\`\`
Do not include explanatory prose outside the code block unless specifically asked.`;

async function callGPT4o(prompt) {
  const { default: OpenAI } = await import("openai");
  const client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  const response = await client.chat.completions.create({
    model: "gpt-4o",
    temperature: 0,
    messages: [
      { role: "system", content: SYSTEM_PROMPT },
      { role: "user",   content: prompt },
    ],
    max_tokens: 2048,
  });
  return response.choices[0].message.content;
}

async function callGemini(prompt) {
  const { GoogleGenerativeAI } = await import("@google/generative-ai");
  const genAI = new GoogleGenerativeAI(process.env.GOOGLE_API_KEY);
  const model = genAI.getGenerativeModel({
    model: "gemini-2.0-flash",
    systemInstruction: SYSTEM_PROMPT,
    generationConfig: { temperature: 0, maxOutputTokens: 2048 },
  });
  const result = await model.generateContent(prompt);
  return result.response.text();
}

try {
  let result;
  if (model === "gpt-4o") {
    result = await callGPT4o(prompt);
  } else if (model === "gemini-2.0-flash") {
    result = await callGemini(prompt);
  } else {
    throw new Error(`Unknown baseline model: ${model}`);
  }
  process.stdout.write(JSON.stringify({ result }) + "\n");
} catch (err) {
  process.stderr.write(`[baseline_shim:${model}] ${err.stack}\n`);
  process.stdout.write(JSON.stringify({ result: null, error: err.message }) + "\n");
  process.exit(1);
}

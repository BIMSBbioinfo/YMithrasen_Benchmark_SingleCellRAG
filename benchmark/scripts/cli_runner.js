/**
 * cli_runner.js
 *
 * Runs a coding agent CLI (Gemini, Claude Code, Codex) non-interactively
 * with a prompt, optionally injecting the recipe as context.
 *
 * Usage: node cli_runner.js <base64-encoded-json>
 * Input JSON: {
 *   agent:       "gemini" | "claude" | "codex",
 *   prompt:      string,
 *   recipe_path: string | null,   // path to recipe.md, or null for vanilla run
 *   work_dir:    string           // temp dir where R code will be executed
 * }
 * Output JSON: { result: string|null, error: string|null }
 */

import { execFile }  from "child_process";
import { promisify } from "util";
import { Buffer }    from "buffer";
import fs            from "fs";
import path          from "path";
import "dotenv/config";

const execFileAsync = promisify(execFile);

// ── helpers ──────────────────────────────────────────────────────────────────

function extractRCode(text) {
  const match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/);
  return match ? match[1].trim() : null;
}

function buildSystemPrompt(recipeContent) {
  const base = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text, comments outside the code block, introductions, or conclusions.
Return ONLY the code block itself.
The files needed for the task are in the current working directory. Use Read10X(data.dir = ".") to load 10x data.`;

  if (!recipeContent) return base;

  return `${base}

---
WORKFLOW RECIPE (use this as your guide):

${recipeContent}
---`;
}

// ── agent callers ─────────────────────────────────────────────────────────────

async function runGeminiCLI(prompt, systemPrompt, workDir) {
  // gemini -p "<prompt>" --system-prompt "<system>" --no-interactive
  // Gemini CLI: https://github.com/google-gemini/gemini-cli
  const fullPrompt = `${systemPrompt}\n\nTask:\n${prompt}`;

  const { stdout } = await execFileAsync(
    "gemini",
    ["--model", "gemini-2.0-flash", "-p", fullPrompt, "--no-interactive"],
    { timeout: 120000, cwd: workDir, env: {
      ...process.env,
      GOOGLE_API_KEY: process.env.GOOGLE_API_KEY
    }}
  );
  return stdout;
}

async function runClaudeCLI(prompt, systemPrompt, workDir) {
  // claude -p "<prompt>" --system "<system>" (Claude Code CLI)
  // Claude Code: https://docs.anthropic.com/en/docs/claude-code
  const { stdout } = await execFileAsync(
    "claude",
    ["-p", prompt, "--system", systemPrompt, "--output-format", "text"],
    { timeout: 120000, cwd: workDir, env: {
      ...process.env,
      ANTHROPIC_API_KEY: process.env.ANTHROPIC_API_KEY
    }}
  );
  return stdout;
}

async function runCodexCLI(prompt, systemPrompt, workDir) {
  // codex "<prompt>" --system "<system>" (OpenAI Codex CLI)
  // Codex CLI: https://github.com/openai/codex
  const fullPrompt = `${systemPrompt}\n\nTask:\n${prompt}`;

  const { stdout } = await execFileAsync(
    "codex",
    ["--full-auto", "-q", fullPrompt],
    { timeout: 120000, cwd: workDir, env: {
      ...process.env,
      OPENAI_API_KEY: process.env.OPENAI_API_KEY
    }}
  );
  return stdout;
}

// ── main ──────────────────────────────────────────────────────────────────────

const payload = JSON.parse(Buffer.from(process.argv[2], "base64").toString("utf8"));
const { agent, prompt, recipe_path, work_dir } = payload;

let result = null;
let error  = null;

try {
  // Load recipe content if provided
  let recipeContent = null;
  if (recipe_path && fs.existsSync(recipe_path)) {
    recipeContent = fs.readFileSync(recipe_path, "utf8");
  }

  const systemPrompt = buildSystemPrompt(recipeContent);

  let raw;
  if      (agent === "gemini") raw = await runGeminiCLI(prompt, systemPrompt, work_dir);
  else if (agent === "claude") raw = await runClaudeCLI(prompt, systemPrompt, work_dir);
  else if (agent === "codex")  raw = await runCodexCLI(prompt, systemPrompt, work_dir);
  else throw new Error(`Unknown agent: ${agent}`);

  result = extractRCode(raw);
  if (!result) throw new Error(`No R code block found in CLI output`);

} catch (e) {
  error = e.message;
  process.stderr.write(`[cli_runner:${agent}] Error: ${error}\n`);
}

process.stdout.write(JSON.stringify({ result, error }) + "\n");

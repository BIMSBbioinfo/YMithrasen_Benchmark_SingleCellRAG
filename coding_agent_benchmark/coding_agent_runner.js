/**
 * coding_agent_runner.js
 *
 * Automated benchmark harness for coding agents on scRNA-seq tasks.
 * Inspired by HumanEval (Pass@k), SWE-bench (repo-grounded eval),
 * Terminal-Bench (end-to-end workflow), and scBench (deterministic grading).
 *
 * NOTE on Pass@k: computePassAtK() below reports whether ANY of the k
 * attempts passed (empirical "success within k attempts"), not HumanEval's
 * formal unbiased Pass@k estimator (which corrects for n>k samples via a
 * combinatorial formula). Report it in the thesis as empirical Pass@k /
 * success-within-k, not as the HumanEval estimator — the two are not
 * numerically interchangeable even though they share a name.
 *
 * Usage:
 *   node coding_agent_runner.js --agent claude
 *   node coding_agent_runner.js --agent codex
 *   node coding_agent_runner.js --agent gemini
 *   node coding_agent_runner.js --agent claude --prompts p01,p02,p03
 *   node coding_agent_runner.js --agent all
 *
 * RESUME: interrupted runs (e.g. Claude Code session/token limits during a
 * long k=3 run) can be resumed with --resume — completed attempts are
 * loaded from results/checkpoints/ and the agent is not called again for
 * them. Resume is OFF by default (a plain invocation always re-runs
 * everything); checkpoints are still written as the run progresses either
 * way, so --resume can pick up from any prior invocation, complete or
 * interrupted. Model outputs are never cached/reused ACROSS attempts (k
 * independence is preserved); only an already-completed, already-graded
 * attempt for the exact same agent+condition+prompt_id+attempt can be
 * skipped — and infrastructure failures (agent_cli_error — no usable
 * response was ever obtained) are never checkpointed, so a session/rate-
 * limit interruption is always retried on --resume rather than treated as
 * done. Checkpoints are automatically invalidated (ignored, re-run) if the
 * prompts file, graders.R, any file under recipes/, or the dataset path
 * changes — see computeBenchmarkVersion().
 */

import "dotenv/config";
import { exec } from "child_process";
import fs   from "fs";
import path from "path";
import crypto from "crypto";
import { fileURLToPath } from "url";

const __dirname = (() => {
  try {
    return path.dirname(fileURLToPath(import.meta.url));
  } catch {
    return process.cwd();
  }
})();

// ── Configuration ──────────────────────────────────────────────────────────────

const CONFIG = {
  // Paths
  repo_root    : path.resolve(__dirname, ".."),
  recipes_dir  : path.resolve(__dirname, "..", "recipes"),
  recipe_file: path.join(__dirname, "..", "recipes", "recipe.md"),
  workflows_dir: path.resolve(__dirname, "..", "recipes"),
  dataset_dir  : process.env.BENCH_DATASET_DIR ||
                 "D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19",
  prompts_file : path.join(__dirname, "prompts", "coding_agent_prompts.jsonl"),
  grader_script: path.join(__dirname, "graders", "graders.R"),
  results_dir  : path.join(__dirname, "results"),
  work_base    : path.join(__dirname, "workdirs"),
  checkpoints_dir: path.join(__dirname, "results", "checkpoints"),
  r_exec       : process.env.R_EXEC || "C:/Program Files/R/R-4.6.1/bin/Rscript.exe",

  // Resume/checkpoint behaviour: OFF by default. A checkpoint is only
  // consulted (and an attempt skipped) when --resume is explicitly passed —
  // safer than resume-by-default, since months later a "new experiment"
  // invocation shouldn't silently reuse old checkpoints. Checkpoints are
  // still WRITTEN after every genuine completed attempt regardless of this
  // flag (see shouldCheckpoint below), so a later --resume run can always
  // pick up from wherever an earlier run left off.
  resume_run   : process.argv.includes("--resume"),

  // Benchmark settings
  k                     : Number(process.env.BENCH_K || 1),
  timeout_ms            : Number(process.env.BENCH_TIMEOUT_MS    || 600000),
  r_timeout_ms          : Number(process.env.BENCH_R_TIMEOUT_MS  || 600000),
  no_recipe             : process.argv.includes("--no-recipe") || process.env.BENCH_NO_RECIPE === "1",
  // Recipe delivery for the non-baseline condition is always progressive
  // disclosure (teaser only, agent discovers ./recipes itself via its own
  // filesystem tools — see buildAgentTeaser). Full-recipe and pointer-card
  // injection modes were removed; recipe_benchmark_runner.js is the frozen
  // full-injection ablation if that comparison is needed again.
  // Context budgeting
  recipe_warning_chars  : Number(process.env.BENCH_RECIPE_WARNING_CHARS    || 16000),
  recipe_max_chars      : Number(process.env.BENCH_RECIPE_MAX_CHARS         || 20000),
  context_window_tokens : Number(process.env.BENCH_CONTEXT_WINDOW_TOKENS   || 128000),
  reserved_output_tokens: Number(process.env.BENCH_MAX_TOKENS               || 12000),
  context_warning_ratio : Number(process.env.BENCH_CONTEXT_WARNING_RATIO   || 0.8),

  // Agent CLI commands (sandbox: recipes only + dataset read-only)
  agents: {
    claude: {
      cmd: (prompt, work_dir) =>
        `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | claude --print"`,
      extract: extractRCode,
    },
    codex: {
      cmd: () => {
        const model = process.env.CODEX_MODEL || "gpt-5.6-terra";
        return (
          `powershell -NoProfile -Command ` +
          `"Get-Content -Raw final_prompt.txt | ` +
          `codex exec ` +
          `-m ${model} ` +
          `--sandbox read-only ` +
          `--skip-git-repo-check ` +
          `--ephemeral ` +
          `--color never ` +
          `-"`
        );
      },
      extract: extractRCode,
    },
    gemini: {
      cmd: () => {
        const model =
          process.env.GEMINI_MODEL ||
          "gemini-3.5-flash";

        // NOTE: -p "." is required to force non-interactive/headless mode.
        // Newer Gemini CLI versions default to interactive mode even when
        // piped via stdin with no positional argument — confirmed by direct
        // testing (stdin-only invocation hung waiting for interactive
        // input instead of exiting). Per `gemini --help`, -p/--prompt's
        // given value is appended with stdin content, so passing a trivial
        // placeholder here is enough to flip the mode switch while the
        // actual prompt content still flows entirely through the existing
        // stdin pipe (not as a command-line argument), preserving the fix
        // for the Windows ~8191-char argument-length truncation bug.
        return (
          `powershell -NoProfile -Command ` +
          `"Get-Content -Raw final_prompt.txt | ` +
          `gemini -m ${model} --output-format text -p '.'"`
        );
      },
      extract: extractRCode,
    },
  },
};

// Single source of truth for the experimental condition label — derived
// directly from CONFIG.no_recipe so it can never drift from what
// buildPromptWithContext() / setupWorkDir() actually did. Saved into both
// per-attempt context_metadata.json and the final summary JSON (see #7).
// Named to match skills_benchmark_runner.js's convention: skill relevance
// AND its primary instruction file (recipe.md) have already been decided
// (oracle-selected) — the agent autonomously selects only which supporting
// files, if any, it needs. This is not a claim about discovering native
// Claude/Codex/Gemini Skills activation, nor about discovering which file
// in ./recipes is authoritative.
function getCondition() {
  return CONFIG.no_recipe ? "no_recipe" : "oracle_activated_progressive_disclosure";
}

// ── Resume / checkpointing ───────────────────────────────────────────────────
// Solves interrupted long k=3 runs (e.g. Claude Code session/token limits)
// without compromising Pass@k: only a COMPLETED, already-graded attempt for
// the exact same agent+condition+prompt_id+attempt can be skipped on resume.
// Model outputs are never cached across attempts — every attempt that
// actually needs to run still calls the agent independently.

// Recursively collects (relative_path, content) for every text file under a
// directory, sorted for a deterministic hash regardless of OS listing order.
// Binary/unreadable files are skipped rather than erroring.
function collectDirTextContent(dir) {
  const out = [];
  function walk(d, rel) {
    let entries;
    try { entries = fs.readdirSync(d, { withFileTypes: true }); } catch { return; }
    for (const e of entries) {
      const full = path.join(d, e.name);
      const relPath = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory()) walk(full, relPath);
      else {
        try { out.push([relPath, fs.readFileSync(full, "utf8")]); } catch { /* skip binary/unreadable */ }
      }
    }
  }
  walk(dir, "");
  out.sort((a, b) => a[0].localeCompare(b[0]));
  return out;
}

let _benchmark_version = null;
function computeBenchmarkVersion() {
  if (_benchmark_version) return _benchmark_version;
  const parts = [CONFIG.prompts_file, CONFIG.grader_script]
    .map(p => { try { return fs.readFileSync(p, "utf8"); } catch { return ""; } });
  // Whole recipes directory (recipe.md + every supporting file under
  // recipes/code/, recipes/scripts/, etc.) — not just recipe.md — so a
  // change to any file the agent could actually read invalidates stale
  // checkpoints.
  for (const [relPath, content] of collectDirTextContent(CONFIG.recipes_dir)) {
    parts.push(`${relPath}\u0000${content}`);
  }
  // Dataset identifier — the path only, not its (large) contents. Switching
  // to a different dataset should invalidate checkpoints even if nothing
  // else changed.
  parts.push(CONFIG.dataset_dir);

  _benchmark_version = crypto.createHash("sha256").update(parts.join("\u0000")).digest("hex").slice(0, 12);
  return _benchmark_version;
}

function checkpointPath(agent_name, prompt_id, attempt) {
  return path.join(CONFIG.checkpoints_dir, agent_name, getCondition(), `${prompt_id}_attempt${attempt}.json`);
}

// Only a genuine completed model attempt is checkpoint-worthy. An
// agent_cli_error means the CLI call itself failed (session limit, rate
// limit, transient error) BEFORE any usable response was obtained — that is
// not a completed attempt, and checkpointing it would make resume skip the
// exact thing it's supposed to retry. grader_fail / r_runtime_error /
// no_code_extracted / timeout all mean the agent DID produce a response
// (which was then graded, or failed to run/parse) — those are genuine
// completed attempts and remain checkpointable.
function shouldCheckpoint(result) {
  return result.failure_type !== "agent_cli_error";
}

function loadCheckpoint(agent_name, prompt_id, attempt) {
  if (!CONFIG.resume_run) return null;
  const p = checkpointPath(agent_name, prompt_id, attempt);
  if (!fs.existsSync(p)) return null;
  try {
    const data = JSON.parse(fs.readFileSync(p, "utf8"));
    if (data.benchmark_version !== computeBenchmarkVersion()) return null; // stale — prompts/recipe/grader/dataset changed
    if (!data.result) return null;
    return data.result;
  } catch {
    return null; // corrupt/partial checkpoint — treat as absent, will be re-run and overwritten
  }
}

function saveCheckpoint(agent_name, prompt_id, attempt, result) {
  const p = checkpointPath(agent_name, prompt_id, attempt);
  fs.mkdirSync(path.dirname(p), { recursive: true });
  fs.writeFileSync(p, JSON.stringify({
    benchmark_version: computeBenchmarkVersion(),
    saved_at: new Date().toISOString(),
    result,
  }, null, 2));
}

// ── Helpers ────────────────────────────────────────────────────────────────────

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);

function log(msg) {
  const ts = new Date().toISOString();
  process.stdout.write(`[${ts}] ${msg}\n`);
}

function compactMarkdown(text) {
  return String(text || "")
    .replace(/\r\n/g, "\n")
    .replace(/[ \t]+$/gm, "")
    .replace(/\n{3,}/g, "\n\n")
    .trim();
}

function cleanPromptWhitespace(text) {
  return String(text || "")
    .replace(/\r\n/g, "\n")
    .replace(/[ \t]+$/gm, "")
    .replace(/\n{3,}/g, "\n\n")
    .trim();
}

function estimateTokens(text) {
  return Math.ceil(String(text || "").length / 3.5);
}

function textStats(text) {
  const value = String(text || "");
  return {
    chars           : value.length,
    bytes           : Buffer.byteLength(value, "utf8"),
    estimated_tokens: estimateTokens(value),
    lines           : value.length === 0 ? 0 : value.split(/\r?\n/).length,
  };
}

function loadJSONL(filepath) {
  log(`Loading prompts from: ${filepath}`);
  if (!fs.existsSync(filepath)) {
    log(`ERROR: prompts file not found at ${filepath}`);
    return [];
  }
  const lines = fs.readFileSync(filepath, "utf8")
    .split("\n").filter(l => l.trim());
  log(`Found ${lines.length} lines in prompts file`);
  const parsed = lines
    .map(l => { try { return JSON.parse(l); } catch(e) { log(`JSON parse error: ${e.message}`); return null; } })
    .filter(Boolean);
  log(`Successfully parsed ${parsed.length} prompts`);
  return parsed;
}

function looksLikeR(code) {
  if (!code) return false;
  return /\b(library|require|Read10X|Read10X_h5|CreateSeuratObject|PercentageFeatureSet|NormalizeData|SCTransform|FindVariableFeatures|ScaleData|RunPCA|RunUMAP|FindNeighbors|FindClusters|FindAllMarkers|FindMarkers|VlnPlot|DimPlot|FeaturePlot|write\.csv|saveRDS)\s*\(/.test(code)
      || /<-\s*/.test(code);
}

function extractRCode(text) {
  if (!text) return null;

  text = String(text)
    .replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "")
    .replace(/<thinking>[\s\S]*?<\/thinking>/gi, "")
    .trim();

  const fenceRegex = /```\s*(?:r|R|\{r[^}]*\})?\s*\r?\n([\s\S]*?)```/g;
  const blocks = [];
  let match;

  while ((match = fenceRegex.exec(text)) !== null) {
    const code = match[1].trim();
    if (looksLikeR(code)) blocks.push(code);
  }

  if (blocks.length > 0) return blocks.join("\n\n");

  if (looksLikeR(text)) return text;

  return null;
}

function loadRecipeDocument() {
  if (!fs.existsSync(CONFIG.recipe_file)) {
    throw new Error(`Recipe file not found: ${CONFIG.recipe_file}`);
  }
  const raw     = fs.readFileSync(CONFIG.recipe_file, "utf8");
  const compact = compactMarkdown(raw);
  const stats   = textStats(compact);

  log(`  [recipe] recipe.md raw_chars=${raw.length} compact_chars=${stats.chars} estimated_tokens=${stats.estimated_tokens}`);

  if (stats.chars > CONFIG.recipe_warning_chars) {
    log(`  [recipe] WARNING: recipe.md is ${stats.chars} chars, above warning threshold ${CONFIG.recipe_warning_chars}`);
  }
  if (stats.chars > CONFIG.recipe_max_chars) {
    throw new Error(`recipe.md is ${stats.chars} compacted chars, above BENCH_RECIPE_MAX_CHARS=${CONFIG.recipe_max_chars}. Shorten the recipe rather than silently truncating it.`);
  }
  return { raw, compact, stats };
}

// Pulls `name:` out of recipe.md's YAML frontmatter (between the first two
// `---` lines) without a full YAML parser — same convention as
// skills_benchmark_runner.js's extractFrontmatterField.
function extractFrontmatterField(raw, field) {
  const frontmatterMatch = raw.match(/^---\n([\s\S]*?)\n---/);
  if (!frontmatterMatch) return null;
  const line = frontmatterMatch[1].split("\n").find(l => l.trim().startsWith(`${field}:`));
  return line ? line.split(":").slice(1).join(":").trim() : null;
}

// Pulls the body paragraph under a given "## Heading" (up to the next "##").
// Same convention as skills_benchmark_runner.js's extractSection.
function extractNamedSection(raw, heading) {
  const re = new RegExp(`##\\s+${heading.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}\\s*\\n([\\s\\S]*?)(?=\\n##\\s|$)`, "i");
  const match = raw.match(re);
  return match ? match[1].trim() : null;
}

// Builds the SHORT metadata-only teaser injected into the prompt for the
// "progressive" condition — mirrors skills_benchmark_runner.js's
// buildSkillTeaser(): the agent gets a skill name + a one-paragraph
// description of what it does, NOT the recipe body or any script content.
//
// The relevant skill AND its primary instruction file (recipe.md) are
// provided by the benchmark condition — the coding agent autonomously
// selects only which, if any, additional supporting files it needs for the
// specific task. This is deliberate: the condition is named
// oracle_activated_progressive_disclosure because skill *selection* is
// already decided; this runner is not testing whether the agent can
// discover which file in ./recipes is authoritative.
//
// The key difference from the plain-API skills condition: these are coding
// agents (Claude Code / Codex / Gemini CLI) with real filesystem access, and
// setupWorkDir() already copies the full ./recipes folder into the agent's
// sandboxed work_dir before it runs. So there's no need for a fake view()
// tool like skills_shim.js provides to the non-agentic API models — the
// agent can just use its own native file-reading tools on ./recipes.
function buildAgentTeaser() {
  const recipe = loadRecipeDocument();
  const name = extractFrontmatterField(recipe.raw, "name") || "a domain-specific single-cell analysis skill";
  const whatItDoes = extractNamedSection(recipe.raw, "What this recipe does") || "(description unavailable)";

  const teaser = [
    "A skill is available and has already been determined to be relevant to the task below.",
    "",
    `Name: ${name}`,
    "",
    `What it does: ${whatItDoes}`,
    "",
    "First read ./recipes/recipe.md in full. This is the primary instruction file for the skill.",
    "Only if recipe.md does not provide enough information for the SPECIFIC task below, read the smallest number of relevant supporting files (e.g. under ./recipes/code/) needed to complete it.",
    "Do not read supporting files unless they are necessary, and do not read every referenced pipeline script.",
    "Once you have enough information, write the final executable R code.",
  ].join("\n");

  return { teaser, recipe };
}

function buildPromptWithContext(prompt_obj) {
  // No-recipe condition: bare task with no guidance
  if (CONFIG.no_recipe) {
    const promptText = cleanPromptWhitespace(`
You are an expert bioinformatician working with single-cell RNA-seq data in R using Seurat.

DATASET PATH:
${CONFIG.dataset_dir}

TASK:
${prompt_obj.prompt}

OUTPUT REQUIREMENTS:
- Output ONLY one fenced R code block starting with \`\`\`r and ending with \`\`\`.
- Do not include explanation before or after the code.
- Do not ask for permission.
- Do not install packages inside the generated code.
- Use set.seed(42).
- Save required outputs to the current working directory.
- Generate the smallest complete executable script.
`);
    log(`  [recipe] condition=${getCondition()}`);
    return {
      promptText,
      contextMetadata: {
        condition               : getCondition(),
        final_prompt            : textStats(promptText),
        configured_context_window_tokens: CONFIG.context_window_tokens,
        reserved_output_tokens  : CONFIG.reserved_output_tokens,
      },
    };
  }

  // Progressive-disclosure condition: the agent gets a short teaser only.
  log(`  [recipe] condition=${getCondition()}`);

  const { teaser, recipe } = buildAgentTeaser();

  const promptText = cleanPromptWhitespace(`
You are an expert bioinformatician working with single-cell RNA-seq data in R using Seurat.

${teaser}

DATASET PATH:
${CONFIG.dataset_dir}

TASK:
${prompt_obj.prompt}

OUTPUT REQUIREMENTS:
- Output ONLY one fenced R code block starting with \`\`\`r and ending with \`\`\`.
- Do not include explanation before or after the code.
- Do not ask for permission.
- Do not install packages inside the generated code.
- Use set.seed(42).
- Save required outputs to the current working directory.
- Generate the smallest complete executable script.
`);

  const promptStats      = textStats(promptText);
  const usableTokens     = Math.max(0, CONFIG.context_window_tokens - CONFIG.reserved_output_tokens);
  const contextRatio     = usableTokens > 0 ? promptStats.estimated_tokens / usableTokens : 1;

  log(`  [context] estimated_tokens=${promptStats.estimated_tokens} budget=${usableTokens} ratio=${(contextRatio * 100).toFixed(1)}%`);

  if (contextRatio >= CONFIG.context_warning_ratio) {
    log(`  [context] WARNING: estimated input uses ${(contextRatio * 100).toFixed(1)}% of token budget`);
  }
  if (promptStats.estimated_tokens > usableTokens) {
    throw new Error(`Estimated prompt size ${promptStats.estimated_tokens} tokens exceeds input budget ${usableTokens}.`);
  }

  return {
    promptText,
    contextMetadata: {
      condition                        : getCondition(),
      teaser                            : textStats(teaser),
      recipe_provenance                 : recipe.stats,
      final_prompt                     : promptStats,
      configured_context_window_tokens : CONFIG.context_window_tokens,
      reserved_output_tokens           : CONFIG.reserved_output_tokens,
      estimated_input_budget_tokens    : usableTokens,
      estimated_context_usage_ratio    : contextRatio,
    },
  };
}

// ── Working directory setup ────────────────────────────────────────────────────

function setupWorkDir(agent, prompt_id, attempt) {
  const work_dir = path.join(
    CONFIG.work_base,
    agent,
    getCondition(),
    `${prompt_id}_attempt${attempt}`
  );

  fs.mkdirSync(work_dir, {
    recursive: true,
  });

  // Remove all outputs from an earlier run.
  cleanWorkDir(work_dir);

  // Copy a fresh recipe directory — ONLY for the progressive-disclosure
  // condition. For --no-recipe, the folder must not exist at all: these are
  // coding agents with general filesystem access, so even an unmentioned
  // ./recipes folder sitting in the work_dir could be discovered and read on
  // its own initiative, silently contaminating the vanilla baseline.
  if (!CONFIG.no_recipe) {
    const recipes_dest = path.join(
      work_dir,
      "recipes"
    );

    fs.cpSync(
      CONFIG.recipes_dir,
      recipes_dest,
      {
        recursive: true,
        force: true,
      }
    );
  }

  return work_dir;
}

function cleanWorkDir(work_dir) {
  const files = fs.readdirSync(work_dir);

  for (const file of files) {
    if (file !== "recipes") {
      fs.rmSync(
        path.join(work_dir, file),
        {
          recursive: true,
          force: true,
        }
      );
    }
  }

  const recipes_dir = path.join(
    work_dir,
    "recipes"
  );

  if (fs.existsSync(recipes_dir)) {
    fs.rmSync(
      recipes_dir,
      {
        recursive: true,
        force: true,
      }
    );
  }
}

// ── OpenRouter API call ────────────────────────────────────────────────────────

async function callOpenRouter(model, prompt_text) {
  const apiKey = process.env.OPENROUTER_API_KEY;
  if (!apiKey) return { stdout: "", stderr: "", error: "OPENROUTER_API_KEY not set" };

  try {
    const response = await fetch("https://openrouter.ai/api/v1/chat/completions", {
      method: "POST",
      headers: {
        "Content-Type":  "application/json",
        "Authorization": `Bearer ${apiKey}`,
        "HTTP-Referer":  "https://github.com/BIMSBbioinfo/YMithrasen_Benchmark_SingleCellRAG",
      },
      body: JSON.stringify({
        model,
        max_tokens: 8192,
        messages: [{ role: "user", content: prompt_text }],
      }),
    });

    if (!response.ok) {
      const body = await response.text();
      return { stdout: "", stderr: body.slice(0, 300), error: `HTTP ${response.status}: ${body.slice(0, 200)}` };
    }

    const data    = await response.json();
    const content = data?.choices?.[0]?.message?.content;
    const text    = Array.isArray(content)
      ? content.map(p => (typeof p === "string" ? p : p?.text || "")).join("")
      : String(content || "");
    return {
      stdout      : text,
      stderr      : "",
      error       : text ? null : "OpenRouter returned empty content",
      model       : data.model || model,
      usage       : data.usage || null,
      finish_reason: data.choices?.[0]?.finish_reason || null,
    };
  } catch (e) {
    return { stdout: "", stderr: e.message, error: e.message };
  }
}

// ── Agent call ─────────────────────────────────────────────────────────────────

async function callAgent(agent_name, prompt_text, work_dir) {
  const agent = CONFIG.agents[agent_name];

  fs.writeFileSync(path.join(work_dir, "final_prompt.txt"), prompt_text, "utf8");
  log(`  [prompt] chars=${prompt_text.length}`);
  log(`  [${agent_name}] calling agent...`);

  // API-mode agents (e.g. gemini via OpenRouter) — no CLI subprocess
  if (agent.api) {
    const result = await callOpenRouter(agent.model, prompt_text);
    fs.writeFileSync(path.join(work_dir, "agent_stdout.txt"), result.stdout || "", "utf8");
    fs.writeFileSync(path.join(work_dir, "agent_stderr.txt"), result.stderr || "", "utf8");
    if (result.error) log(`  [${agent_name}] api error: ${result.error.slice(0, 180)}`);
    return result;
  }

  // CLI-mode agents (claude, codex)
  return new Promise((resolve) => {
    const cmd = agent.cmd(prompt_text, work_dir);
    const env = { ...process.env };

    exec(cmd, {
      timeout  : CONFIG.timeout_ms,
      cwd      : work_dir,
      env,
      maxBuffer: 20 * 1024 * 1024,
    }, (error, stdout, stderr) => {
      fs.writeFileSync(path.join(work_dir, "agent_stdout.txt"), stdout || "", "utf8");
      fs.writeFileSync(path.join(work_dir, "agent_stderr.txt"), stderr || "", "utf8");
      if (error) log(`  [${agent_name}] agent error: ${error.message?.slice(0, 180)}`);
      resolve({ stdout: stdout || "", stderr: stderr || "", error: error?.message || null });
    });
  });
}

// ── R execution ───────────────────────────────────────────────────────────────

function executeRCode(r_code, work_dir, prompt_id) {
  return new Promise((resolve) => {
    const r_file = path.join(work_dir, `${prompt_id}_generated.R`);
    fs.writeFileSync(r_file, r_code, "utf8");

    log(`  [R] executing ${path.basename(r_file)}`);

    exec(
      `"${CONFIG.r_exec}" "${r_file}" 2>&1`,
      { timeout: CONFIG.r_timeout_ms, cwd: work_dir, maxBuffer: 5 * 1024 * 1024 },
      (error, stdout) => {
        const success = !error || error.code === 0;
        resolve({
          success,
          stdout: stdout || "",
          error : error?.message || null,
          r_file,
        });
      }
    );
  });
}

// ── Grader ─────────────────────────────────────────────────────────────────────

function runGrader(prompt_id, work_dir) {
  return new Promise((resolve) => {
    const cmd = `"${CONFIG.r_exec}" "${CONFIG.grader_script}" "${prompt_id}" "${work_dir}"`;
    exec(cmd, { timeout: 30000, maxBuffer: 1024 * 1024 }, (error, stdout) => {
      try {
        const lines  = (stdout || "").trim().split("\n").filter(l => l.startsWith("{"));
        const result = JSON.parse(lines.pop() || "{}");
        resolve(result);
      } catch {
        resolve({
          prompt_id,
          pass           : false,
          score          : 0.0,
          reason         : "Grader parse error: " + (error?.message || "unknown"),
        });
      }
    });
  });
}

// ── Single prompt evaluation ───────────────────────────────────────────────────
// NOTE on `usage`: the claude/codex CLI agents (see CONFIG.agents) are invoked
// via plain stdout/stderr subprocess calls, not a structured API response, so
// their `agent_result.usage` (and therefore every result object's `usage`
// field below) will generally be null — token counts and any provider-side
// prompt-caching effects are NOT captured for those two agents. Only the
// OpenRouter-API-mode path (callOpenRouter) returns real `usage`. If cost /
// caching analysis becomes an explicit thesis result, this runner does not
// currently support it for CLI agents; otherwise it's enough to state that
// each attempt was an independent invocation and caching was not controlled
// or measured as a variable.
async function evaluatePrompt(agent_name, prompt_obj, attempt = 1) {
  const { id: prompt_id } = prompt_obj;
  log(`\n  [${agent_name}/${prompt_id}] attempt ${attempt}`);

  const work_dir = setupWorkDir(agent_name, prompt_id, attempt);
  const { promptText: prompt_text, contextMetadata } = buildPromptWithContext(prompt_obj);

  fs.writeFileSync(path.join(work_dir, "context_metadata.json"), JSON.stringify(contextMetadata, null, 2), "utf8");

  const t_start = Date.now();

  // 1. Call agent
  const agent_result = await callAgent(agent_name, prompt_text, work_dir);

  const actualPromptTokens = agent_result.usage?.prompt_tokens ?? null;
  if (actualPromptTokens !== null) {
    const estimated = contextMetadata?.final_prompt?.estimated_tokens || 0;
    log(`  [context] estimated_input_tokens=${estimated} actual_input_tokens=${actualPromptTokens}`);
    if (estimated > 0 && Math.abs(actualPromptTokens - estimated) / estimated > 0.25) {
      log(`  [context] NOTE: token estimate differs from provider count by more than 25%`);
    }
  }
  const t_agent      = Date.now() - t_start;
  log(`  [${agent_name}/${prompt_id}] agent responded in ${(t_agent/1000).toFixed(1)}s`);

  // Handle agent CLI errors (no output + error message)
  if (agent_result.error && !agent_result.stdout) {
    return {
      agent: agent_name,
      prompt_id,
      category         : prompt_obj.category || prompt_obj.id?.split('_')[0],
      difficulty: prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : Date.now() - t_start,
      r_executed       : false,
      r_success        : false,
      grader_pass      : false,
      grader_score     : 0.0,
      failure_type     : "agent_cli_error",
      reason           : agent_result.error,
      context          : contextMetadata,
      model            : agent_result.model  || CONFIG.agents[agent_name]?.model || null,
      usage            : agent_result.usage  || null,
      finish_reason    : agent_result.finish_reason || null,
      work_dir,
    };
  }

  // 2. Extract R code
  const r_code = extractRCode(agent_result.stdout);

  if (!r_code) {
    log(`  [${agent_name}/${prompt_id}] no R code extracted`);
    return {
      agent: agent_name,
      prompt_id,
      category         : prompt_obj.category || prompt_obj.id?.split('_')[0],
      difficulty: prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : Date.now() - t_start,
      r_executed       : false,
      r_success        : false,
      grader_pass      : false,
      grader_score     : 0.0,
      failure_type     : "no_code_extracted",
      reason           : "No R code in agent output. Inspect agent_stdout.txt and agent_stderr.txt in the workdir.",
      context          : contextMetadata,
      model            : agent_result.model  || CONFIG.agents[agent_name]?.model || null,
      usage            : agent_result.usage  || null,
      finish_reason    : agent_result.finish_reason || null,
      work_dir,
    };
  }

  fs.writeFileSync(path.join(work_dir, `${prompt_id}_extracted.R`), r_code, "utf8");

  // 3. Execute R code
  const r_result = await executeRCode(r_code, work_dir, prompt_id);
  fs.writeFileSync(path.join(work_dir, `${prompt_id}_r_stdout.txt`), r_result.stdout || "", "utf8");

  const t_total  = Date.now() - t_start;
  log(`  [${agent_name}/${prompt_id}] R ${r_result.success ? "✓" : "✗"} (${(t_total/1000).toFixed(1)}s)`);

  if (!r_result.success) {
    return {
      agent          : agent_name,
      prompt_id,
      category       : prompt_obj.category || prompt_obj.id?.split("_")[0],
      difficulty     : prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : t_total,
      r_executed       : true,
      r_success        : false,
      r_error          : r_result.error,
      grader_pass      : false,
      grader_score     : 0.0,
      failure_type     : r_result.error && /timed out|timeout/i.test(r_result.error) ? "timeout" : "r_runtime_error",
      reason           : r_result.error || "R execution failed. Inspect *_r_stdout.txt in the workdir.",
      context          : contextMetadata,
      model            : agent_result.model  || CONFIG.agents[agent_name]?.model || null,
      usage            : agent_result.usage  || null,
      finish_reason    : agent_result.finish_reason || null,
      work_dir,
    };
  }
  // 4. Run grader
  const grader_result = await runGrader(prompt_id, work_dir);
  log(`  [${agent_name}/${prompt_id}] grader: pass=${grader_result.pass} score=${grader_result.score?.toFixed(2)}`);

  return {
    agent          : agent_name,
    prompt_id,
    category       : prompt_obj.category || prompt_obj.id?.split("_")[0],
    difficulty     : prompt_obj.difficulty,
    attempt,
    agent_latency_ms : t_agent,
    total_ms         : t_total,
    r_executed       : true,
    r_success        : r_result.success,
    r_error          : r_result.error,
    grader_pass      : grader_result.pass   || false,
    grader_score     : grader_result.score  || 0.0,
    failure_type     : grader_result.pass ? "pass" : "grader_fail",
    reason           : grader_result.reason || "",
    context          : contextMetadata,
    model            : agent_result.model  || CONFIG.agents[agent_name]?.model || null,
    usage            : agent_result.usage  || null,
    finish_reason    : agent_result.finish_reason || null,
    work_dir,
  };
}

// ── Pass@k calculation ─────────────────────────────────────────────────────────

// This is empirical "success within k attempts" (did any of the k independent
// attempts pass), NOT HumanEval's unbiased Pass@k estimator. Report it as such.
function computePassAtK(attempts) {
  // Pass@k: 1 if any attempt passes, 0 otherwise
  return attempts.some(a => a.grader_pass) ? 1.0 : 0.0;
}

// ── Main benchmark runner ──────────────────────────────────────────────────────

// Pure computation — no logging, no side effects. Called after every
// attempt (cheap even at k=3 × 30 prompts) so results/<file>.json always
// reflects the latest completed state, and once more at the end for the
// final logged summary.
function computeSummary(agent_name, selected, all_results) {
  const prompt_groups = {};
  for (const r of all_results) {
    if (!prompt_groups[r.prompt_id]) prompt_groups[r.prompt_id] = [];
    prompt_groups[r.prompt_id].push(r);
  }
  const prompt_passes  = Object.values(prompt_groups).map(g => computePassAtK(g));
  const pass_rate      = prompt_passes.reduce((s, v) => s + v, 0) / Math.max(prompt_passes.length, 1);
  const r_success_rate = all_results.filter(r => r.r_success).length  / Math.max(all_results.length, 1);
  const avg_score      = all_results.reduce((s, r) => s + (r.grader_score   || 0), 0) / Math.max(all_results.length, 1);
  const avg_latency    = all_results.reduce((s, r) => s + (r.agent_latency_ms || 0), 0) / Math.max(all_results.length, 1);
  const categories     = [...new Set(all_results.map(r => r.category).filter(Boolean))];

  return {
    run_id, agent: agent_name,
    condition   : getCondition(),
    benchmark_version   : computeBenchmarkVersion(),
    timestamp   : new Date().toISOString(),
    n_prompts   : selected.length,
    k           : CONFIG.k,
    attempts_per_prompt : CONFIG.k,
    n_attempts           : all_results.length,
    complete              : all_results.length >= selected.length * CONFIG.k,
    metrics     : { pass_at_k: pass_rate, r_success_rate, avg_score, avg_latency_ms: avg_latency },
    by_category : Object.fromEntries(categories.map(cat => {
      const cat_groups = Object.values(prompt_groups).filter(g => g[0].category === cat);
      return [cat, cat_groups.length ? cat_groups.reduce((s, g) => s + computePassAtK(g), 0) / cat_groups.length : null];
    })),
    by_difficulty : Object.fromEntries(["simple", "medium", "hard"].map(diff => {
      const diff_groups = Object.values(prompt_groups).filter(g =>
        (g[0].difficulty === diff) || (g[0].prompt_id && g[0].prompt_id.startsWith(diff.slice(0, 3))));
      return [diff, diff_groups.length ? diff_groups.reduce((s, g) => s + computePassAtK(g), 0) / diff_groups.length : null];
    })),
    results: all_results,
  };
}

function logSummary(agent_name, summary) {
  log(`\n${"─".repeat(70)}`);
  log(`BENCHMARK SUMMARY — ${agent_name}`);
  log(`${"─".repeat(70)}`);
  log(`Pass@${CONFIG.k} (strict):     ${(summary.metrics.pass_at_k * 100).toFixed(1)}%`);
  log(`Workflow success (R):  ${(summary.metrics.r_success_rate * 100).toFixed(1)}%`);
  log(`Avg grader score:      ${summary.metrics.avg_score.toFixed(3)}`);
  log(`Avg agent latency:     ${(summary.metrics.avg_latency_ms / 1000).toFixed(1)}s`);

  if (Object.keys(summary.by_category).length > 0) {
    log(`\nBy category (prompt-level Pass@${CONFIG.k}):`);
    for (const [cat, v] of Object.entries(summary.by_category)) {
      if (v === null) continue;
      log(`  ${cat.padEnd(20)} ${(v * 100).toFixed(1)}%`);
    }
  }

  log(`\nBy difficulty (prompt-level Pass@${CONFIG.k}):`);
  for (const [diff, v] of Object.entries(summary.by_difficulty)) {
    if (v === null) continue;
    log(`  ${diff.padEnd(10)} ${(v * 100).toFixed(1)}%`);
  }

  log(`\nBy failure type:`);
  const failure_types = [...new Set(summary.results.map(r => r.failure_type || (r.grader_pass ? "pass" : "unknown")))];
  for (const ft of failure_types) {
    const n = summary.results.filter(r => (r.failure_type || (r.grader_pass ? "pass" : "unknown")) === ft).length;
    log(`  ${ft.padEnd(20)} ${n}`);
  }
}

async function runBenchmark(agent_name, prompt_ids = null) {
  log(`\n${"═".repeat(70)}`);
  log(`Coding Agent Benchmark  run_id=${run_id}`);
  log(`Agent: ${agent_name}  k=${CONFIG.k}  condition=${getCondition()}  ${CONFIG.resume_run ? "(--resume: consulting checkpoints)" : "(fresh run — no checkpoints consulted)"}`);
  log(`${"═".repeat(70)}`);

  const prompts = loadJSONL(CONFIG.prompts_file);
  const selected = prompt_ids
    ? prompts.filter(p => prompt_ids.includes(p.id))
    : prompts;

  log(`Loaded ${selected.length} prompts`);
  fs.mkdirSync(CONFIG.results_dir, { recursive: true });
  fs.mkdirSync(CONFIG.work_base,   { recursive: true });

  const out_file = path.join(CONFIG.results_dir, `coding_agent_${agent_name}_${run_id}.json`);
  const all_results = [];

  for (const prompt_obj of selected) {
    const attempts = [];
    for (let k = 1; k <= CONFIG.k; k++) {
      let result = loadCheckpoint(agent_name, prompt_obj.id, k);
      if (result) {
        log(`  [${agent_name}/${prompt_obj.id}] attempt ${k} — resumed from checkpoint (skipping agent call)`);
      } else {
        result = await evaluatePrompt(agent_name, prompt_obj, k);

        if (!shouldCheckpoint(result)) {
          // Infrastructure failure (e.g. Claude session/rate-limit error) —
          // not a model attempt, not counted toward Pass@k, not checkpointed.
          // Stop the run entirely rather than silently moving to the next
          // prompt/attempt, so this can never be mistaken for a genuine
          // failed benchmark attempt.
          log(`  [checkpoint] infrastructure failure (${result.failure_type}) at ${prompt_obj.id} attempt ${k} — attempt NOT counted, run stopping. Resume with --resume once resolved.`);

          const partial_summary = computeSummary(agent_name, selected, all_results);
          fs.writeFileSync(out_file, JSON.stringify(partial_summary, null, 2));

          throw new Error(`Benchmark interrupted at ${prompt_obj.id} attempt ${k}: ${result.reason}`);
        }

        saveCheckpoint(agent_name, prompt_obj.id, k, result);
      }
      attempts.push(result);
      all_results.push(result);

      // Incremental save after every single attempt — a crash mid-run still
      // leaves a valid, up-to-date results JSON rather than nothing at all.
      const partial_summary = computeSummary(agent_name, selected, all_results);
      fs.writeFileSync(out_file, JSON.stringify(partial_summary, null, 2));
    }

    const pass_at_k = computePassAtK(attempts);
    log(`  [${prompt_obj.id}] Pass@${CONFIG.k}=${pass_at_k.toFixed(3)} (${attempts.length}/${CONFIG.k} attempts run)`);
  }

  const summary = computeSummary(agent_name, selected, all_results);
  logSummary(agent_name, summary);
  fs.writeFileSync(out_file, JSON.stringify(summary, null, 2));
  log(`\nResults saved → ${out_file}`);
  log(`✅ Done.`);

  return summary;
}

// ── CLI entry ──────────────────────────────────────────────────────────────────

// Parse CLI args — works both with `node file.js` and `node -e "import(...)"`
const raw_args = process.argv.slice(2);
// Also check env vars as fallback
const agent_idx   = raw_args.indexOf("--agent");
const agent_arg   = agent_idx >= 0 ? raw_args[agent_idx + 1] : (process.env.BENCH_AGENT || "claude");

const prompts_idx = raw_args.indexOf("--prompts");
const prompts_arg = prompts_idx >= 0 ? raw_args[prompts_idx + 1] : (process.env.BENCH_PROMPTS || null);
const prompt_ids  = prompts_arg ? prompts_arg.split(",").map(p => p.trim()) : null;

if (agent_arg === "all") {
  for (const agent of ["claude", "codex", "gemini"]) {
    await runBenchmark(agent, prompt_ids);
  }
} else {
  await runBenchmark(agent_arg, prompt_ids);
}
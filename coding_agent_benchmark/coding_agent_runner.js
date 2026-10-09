// coding_agent_runner.js
// Benchmarks Claude Code, Codex CLI and Gemini CLI on the scRNA-seq prompts.
// Each attempt: build the prompt, call the agent in a fresh working directory,
// extract the R code, run it, then grade the outputs with graders/graders.R.
//
// Two conditions:
//   default      With Skills: a short teaser in the prompt and a copy of recipes/
//                in the working directory (condition oracle_activated_progressive_disclosure)
//   --no-recipe  Without Skills: the task prompt only (condition no_recipe)
//
// Usage (from the repo root; k comes from BENCH_K, default 1):
//   node coding_agent_benchmark/coding_agent_runner.js --agent claude
//   node coding_agent_benchmark/coding_agent_runner.js --agent codex --no-recipe
//   node coding_agent_benchmark/coding_agent_runner.js --agent all --prompts simple_01,hard_04
//   node coding_agent_benchmark/coding_agent_runner.js --agent claude --resume
//
// Pass@k here means "at least one of the k attempts passed", not the HumanEval estimator.
// Interrupted runs can be continued with --resume: completed attempts are read from
// results/checkpoints/ and CLI failures are always retried.
// Needs the claude, codex and gemini CLIs on the PATH.

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

// ---- Configuration ----

// The prompt is piped to the agent CLI from final_prompt.txt:
// through PowerShell on Windows, through sh on Linux and macOS.
function pipePrompt(cli) {
  return process.platform === "win32"
    ? `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | ${cli}"`
    : `cat final_prompt.txt | ${cli}`;
}

const CONFIG = {
  // Paths
  repo_root    : path.resolve(__dirname, ".."),
  recipes_dir  : path.resolve(__dirname, "..", "recipes"),
  recipe_file: path.join(__dirname, "..", "recipes", "recipe.md"),
  workflows_dir: path.resolve(__dirname, "..", "recipes"),
  dataset_dir  : process.env.BENCH_DATASET_DIR ||
                 path.resolve(__dirname, "..", "data", "pbmc3k", "filtered_gene_bc_matrices", "hg19").replace(/\\/g, "/"),
  prompts_file : path.join(__dirname, "prompts", "coding_agent_prompts.jsonl"),
  grader_script: path.join(__dirname, "graders", "graders.R"),
  results_dir  : path.join(__dirname, "results"),
  work_base    : path.join(__dirname, "workdirs"),
  checkpoints_dir: path.join(__dirname, "results", "checkpoints"),
  r_exec       : process.env.R_EXEC || "Rscript", // Rscript from PATH unless R_EXEC is set

  // checkpoints are always written but only read with --resume
  resume_run   : process.argv.includes("--resume"),

  // Benchmark settings
  k                     : Number(process.env.BENCH_K || 1),
  timeout_ms            : Number(process.env.BENCH_TIMEOUT_MS    || 600000),
  r_timeout_ms          : Number(process.env.BENCH_R_TIMEOUT_MS  || 600000),
  no_recipe             : process.argv.includes("--no-recipe") || process.env.BENCH_NO_RECIPE === "1",
  // the skill is always delivered as teaser + files; see buildAgentTeaser
  recipe_warning_chars  : Number(process.env.BENCH_RECIPE_WARNING_CHARS    || 16000),
  recipe_max_chars      : Number(process.env.BENCH_RECIPE_MAX_CHARS         || 20000),
  context_window_tokens : Number(process.env.BENCH_CONTEXT_WINDOW_TOKENS   || 128000),
  reserved_output_tokens: Number(process.env.BENCH_MAX_TOKENS               || 12000),
  context_warning_ratio : Number(process.env.BENCH_CONTEXT_WARNING_RATIO   || 0.8),

  // agent CLI commands; the prompt is piped in from final_prompt.txt
  agents: {
    claude: {
      cmd: (prompt, work_dir) => pipePrompt("claude --print"),
      extract: extractRCode,
    },
    codex: {
      cmd: () => {
        const model = process.env.CODEX_MODEL || "gpt-5.6-terra";
        return pipePrompt(
          `codex exec ` +
          `-m ${model} ` +
          `--sandbox read-only ` +
          `--skip-git-repo-check ` +
          `--ephemeral ` +
          `--color never ` +
          `-`
        );
      },
      extract: extractRCode,
    },
    gemini: {
      cmd: () => {
        const model =
          process.env.GEMINI_MODEL ||
          "gemini-3.5-flash";

        // -p '.' forces headless mode; the real prompt still comes in on stdin
        return pipePrompt(`gemini -m ${model} --output-format text -p '.'`);
      },
      extract: extractRCode,
    },
  },
};

// condition label saved in every result file
function getCondition() {
  return CONFIG.no_recipe ? "no_recipe" : "oracle_activated_progressive_disclosure";
}

// ---- Resume / checkpointing ----
// Only a completed, graded attempt for the same agent, condition, prompt and attempt
// number is reused. Model output is never shared between attempts.

// [relative path, content] for every text file under dir, sorted by path
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
// Hash of the prompts file, graders.R, everything under recipes/ and the dataset path.
// A checkpoint from a different version is ignored.
function computeBenchmarkVersion() {
  if (_benchmark_version) return _benchmark_version;
  const parts = [CONFIG.prompts_file, CONFIG.grader_script]
    .map(p => { try { return fs.readFileSync(p, "utf8"); } catch { return ""; } });
  // the whole recipes/ folder, since the agent can read any file in it
  for (const [relPath, content] of collectDirTextContent(CONFIG.recipes_dir)) {
    parts.push(`${relPath}\u0000${content}`);
  }
  // dataset path only, not its contents
  parts.push(CONFIG.dataset_dir);

  _benchmark_version = crypto.createHash("sha256").update(parts.join("\u0000")).digest("hex").slice(0, 12);
  return _benchmark_version;
}

function checkpointPath(agent_name, prompt_id, attempt) {
  return path.join(CONFIG.checkpoints_dir, agent_name, getCondition(), `${prompt_id}_attempt${attempt}.json`);
}

// agent_cli_error means the CLI call itself failed (session or rate limit), so there is
// nothing to checkpoint and the attempt must be retried.
function shouldCheckpoint(result) {
  return result.failure_type !== "agent_cli_error";
}

function loadCheckpoint(agent_name, prompt_id, attempt) {
  if (!CONFIG.resume_run) return null;
  const p = checkpointPath(agent_name, prompt_id, attempt);
  if (!fs.existsSync(p)) return null;
  try {
    const data = JSON.parse(fs.readFileSync(p, "utf8"));
    if (data.benchmark_version !== computeBenchmarkVersion()) return null; // stale: prompts, recipe, grader or dataset changed
    if (!data.result) return null;
    return data.result;
  } catch {
    return null; // corrupt checkpoint: treat as missing
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

// ---- Helpers ----

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

// rough estimate, for logging only
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

// Drops <scratchpad>/<thinking> blocks, then joins all fenced blocks that look like R.
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

// reads a field from recipe.md's YAML frontmatter
function extractFrontmatterField(raw, field) {
  const frontmatterMatch = raw.match(/^---\n([\s\S]*?)\n---/);
  if (!frontmatterMatch) return null;
  const line = frontmatterMatch[1].split("\n").find(l => l.trim().startsWith(`${field}:`));
  return line ? line.split(":").slice(1).join(":").trim() : null;
}

// body of a "## Heading" section, up to the next "##"
function extractNamedSection(raw, heading) {
  const re = new RegExp(`##\\s+${heading.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")}\\s*\\n([\\s\\S]*?)(?=\\n##\\s|$)`, "i");
  const match = raw.match(re);
  return match ? match[1].trim() : null;
}

// Teaser for the With Skills condition: skill name, what it does, and where recipe.md is.
// The agent reads the files in ./recipes with its own tools.
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

// Full prompt for one task: role, teaser (With Skills only), dataset path, task, output rules.
function buildPromptWithContext(prompt_obj) {
  // Without Skills: task only
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

  // With Skills: task plus teaser
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

// ---- Working directory setup ----

// Fresh working directory for one attempt.
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

  // recipes/ is copied only in the With Skills condition. Without Skills it must not
  // exist at all, or the agent could find it on its own.
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

// ---- OpenRouter API call ----

// Only for agents declared with api: true. None of the three configured agents use it.
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

// ---- Agent call ----

// Writes the prompt to final_prompt.txt and pipes it to the agent CLI
// (avoids the Windows command-line length limit).
async function callAgent(agent_name, prompt_text, work_dir) {
  const agent = CONFIG.agents[agent_name];

  fs.writeFileSync(path.join(work_dir, "final_prompt.txt"), prompt_text, "utf8");
  log(`  [prompt] chars=${prompt_text.length}`);
  log(`  [${agent_name}] calling agent...`);

  // API-mode agents: no CLI subprocess
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

// ---- R execution ----

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

// ---- Grader ----

// Runs graders.R on the working directory and parses its JSON verdict.
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

// ---- Single prompt evaluation ----
// Token usage is not captured for CLI agents, so `usage` is null in their results.
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

// ---- Pass@k calculation ----

// empirical Pass@k: did any of the k attempts pass
function computePassAtK(attempts) {
  return attempts.some(a => a.grader_pass) ? 1.0 : 0.0;
}

// ---- Main benchmark runner ----

// Builds the result-file object from the attempts so far. Called after every attempt.
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

// Runs k attempts per prompt for one agent and writes the results file.
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
          // CLI failure: not a model attempt, so stop the run instead of counting it. Use --resume to continue.
          log(`  [checkpoint] infrastructure failure (${result.failure_type}) at ${prompt_obj.id} attempt ${k} — attempt NOT counted, run stopping. Resume with --resume once resolved.`);

          const partial_summary = computeSummary(agent_name, selected, all_results);
          fs.writeFileSync(out_file, JSON.stringify(partial_summary, null, 2));

          throw new Error(`Benchmark interrupted at ${prompt_obj.id} attempt ${k}: ${result.reason}`);
        }

        saveCheckpoint(agent_name, prompt_obj.id, k, result);
      }
      attempts.push(result);
      all_results.push(result);

      // save after every attempt so a crash still leaves a valid results file
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

// ---- CLI entry ----

// CLI args, with env vars as fallback
const raw_args = process.argv.slice(2);
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
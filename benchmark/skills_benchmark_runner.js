/**
 * skills_benchmark_runner.js
 *
 * Sibling of recipe_benchmark_runner.js. Where that script pastes the
 * ENTIRE recipe.md text into the system prompt up front (then allows
 * fetching pre-announced supporting files), this script instead SIMULATES
 * the shape of Claude's Skills progressive disclosure, as closely as a bare
 * chat-completion API call reasonably allows:
 *   - Only a SHORT metadata teaser (name + what-the-skill-does, built by
 *     buildSkillTeaser() below) is visible up front — not the full recipe.
 *   - The model must discover and read recipe.md itself via a general
 *     `view(path)` tool (list-a-directory-or-read-a-file, sandboxed to the
 *     recipes/ folder but not restricted to a fixed filename list) before
 *     it has any procedural instructions at all.
 *   - The same tool is then used again to browse/read supporting files
 *     (code/*.R, workflow JSONs) — the model finds these itself rather than
 *     being told a filename in advance.
 * This is a reasonable simulation of the progressive-disclosure PATTERN,
 * not a claim of equivalence to Claude's actual internal skill-loading
 * implementation, which this project has no visibility into. The `view`
 * tool itself and the multi-turn tool loop live in skills_shim.js, not
 * here — this script's job is building the teaser, running the same proven
 * R-sandbox/scoring/preflight machinery as recipe_benchmark_runner.js, and
 * keeping output files clearly separate (skills_results_*.json /
 * skills_comparison_*.json) so this condition never overwrites or gets
 * confused with the other one.
 *
 * IMPORTANT DISCLOSURE — forced final synthesis is NOT documented Claude
 * Skills behavior. Anthropic's tool-use docs describe an ordinary
 * request/respond loop with no arbitrary turn cap or mandatory
 * tools-disabled final turn. The bounded-exploration-then-forced-answer
 * mechanism in skills_shim.js (MAX_TOOL_TURNS, forceFinal*Answer()) is THIS
 * BENCHMARK'S OWN cost/robustness policy, added so a model that's genuinely
 * still making progress (reading new files, not looping) doesn't get
 * discarded as a hard failure purely because of an arbitrary turn cap. The
 * condition label reflects this precisely: "bounded_exploration_
 * forced_synthesis", not a claim about how Claude Skills actually work.
 * `tool_budget_exhausted` is recorded per attempt specifically so a model
 * that needed forced synthesis is distinguishable from one that read a few
 * relevant files and stopped voluntarily — report the `Forced` count
 * alongside Pass@k when interpreting results, not just the headline score.
 *
 * One deliberate simplification kept from recipe_benchmark_runner.js: skill
 * SELECTION itself is still skipped — the teaser tells the model this skill
 * IS the relevant one (oracle-selected), rather than making it choose among
 * several candidate skills. That's a further, separate extension.
 *
 * Condition modeled: "Oracle-selected skill with metadata-only upfront
 * context and on-demand access to the main recipe and supporting files
 * through a sandboxed view tool" — distinct from:
 *   - benchmark_runner.js's vanilla baseline (no skill context, no tools)
 *   - recipe_benchmark_runner.js's condition (full recipe.md pre-pasted,
 *     narrow enum-constrained tool for supporting files only)
 *   - coding_agent_runner.js's CLI-agent "with-recipe" condition (a full
 *     coding agent with general filesystem access, not a scoped tool)
 *
 * recipe.md, once fetched via view("recipe.md"), arrives as RAW, UNMODIFIED
 * TEXT (no whitespace normalisation, no summarisation) — same principle as
 * recipe_benchmark_runner.js: summarising would measure a different
 * condition. viewPath() in skills_shim.js does no text transformation.
 *
 * IMPORTANT — injection channel: the teaser is sent to skills_shim.js as a
 * separate `systemContext` field, NOT concatenated into the user/task
 * prompt string, same convention as recipe_benchmark_runner.js. The task
 * prompt (`prompt`) is sent completely unchanged from the vanilla baseline.
 *
 * Sandbox: identical to benchmark_runner.js / recipe_benchmark_runner.js —
 * each attempt's R code is written to its own temp dir under
 * CONFIG.results_dir and executed via `Rscript --vanilla` as a child
 * process with a hard timeout. Temp dirs are cleaned up on success, left on
 * disk on failure for post-hoc debugging.
 *
 * This script produces output in the SAME JSON schema as benchmark_runner.js
 * / recipe_benchmark_runner.js (run_id / config / results{model:{passK,
 * quality}}), which is what makes a direct diff against the frozen vanilla
 * baseline meaningful — same prompts, same k, same validation logic, same
 * grading code. The only variable being changed is HOW the skill is
 * accessed.
 *
 * Usage:
 *   node skills_benchmark_runner.js
 *   node skills_benchmark_runner.js --model claude-sonnet-5,gpt-5.5
 *   node skills_benchmark_runner.js --k 3 --baseline results/results_2026-07-19T20-16-10-541.json
 *   node skills_benchmark_runner.js --prompts simple_01,medium_04 --allow-partial-prompts
 *   node skills_benchmark_runner.js --skip-comparison
 *
 * RECOMMENDED before a full run: do a single-model, single-prompt dry run
 * and check the preflight log for "teaser confirmed applied" AND "view tool
 * confirmed working (list + read)" — separate mechanisms that can fail
 * independently. Also inspect the saved skills_results_*.json's per-attempt
 * `tool_calls_made` field and the console's "ToolUse" column to see whether
 * the model actually discovered and read recipe.md + supporting files on
 * real prompts, not just the smoke test's explicitly-instructed one:
 *   node skills_benchmark_runner.js --model claude-sonnet-5 --prompts simple_01 --k 1 --skip-comparison
 */

import "dotenv/config";
import fs     from "fs";
import path   from "path";
import crypto from "crypto";
import { execFile, exec, execFileSync } from "child_process";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

function sha256(text) {
  return crypto.createHash("sha256").update(String(text || ""), "utf8").digest("hex");
}

// ── Config ────────────────────────────────────────────────────────────────────

function argVal(flag, fallback) {
  const v = process.argv.find((a, i) => process.argv[i - 1] === flag);
  return v !== undefined ? v : fallback;
}

const CONFIG = {
  k:            parseInt(argVal("--k", "3")),
  model:        argVal("--model", "all"),
  prompts:      argVal("--prompts", "all"),
  allowPartialPrompts: process.argv.includes("--allow-partial-prompts"),
  skipComparison:      process.argv.includes("--skip-comparison"),

  r_exec:       process.env.R_EXEC || "Rscript",
  dataset_dir:  (process.env.BENCH_DATASET_DIR ||
    "C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19"
  ).replace(/\\/g, "/"),

  results_dir:  path.join(__dirname, "results"),
  timeout_ms:   600000, // 10 min, same as benchmark_runner.js

  prompts_file: path.join(__dirname, "prompts/benchmark_prompts.jsonl"),

  // Repo layout: benchmark/recipe_benchmark_runner.js, recipes/recipe.md is a
  // sibling of benchmark/ — same relative path convention as
  // coding_agent_runner.js's CONFIG.recipe_file.
  recipe_file:  path.resolve(__dirname, "..", "recipes", "recipe.md"),
  recipe_warning_chars: Number(process.env.BENCH_RECIPE_WARNING_CHARS || 16000),
  recipe_max_chars:     Number(process.env.BENCH_RECIPE_MAX_CHARS     || 20000),

  // Frozen vanilla-baseline results file this run gets diffed against.
  baseline_file: argVal("--baseline",
    path.join(__dirname, "results", "results_2026-07-19T20-16-10-541.json")),
};

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 23);
fs.mkdirSync(CONFIG.results_dir, { recursive: true });

// ── Logging ───────────────────────────────────────────────────────────────────

function log(msg) {
  console.log(`[${new Date().toISOString()}] ${msg}`);
}

// ── Load prompts ──────────────────────────────────────────────────────────────

function loadJSONL(file) {
  return fs.readFileSync(file, "utf8")
    .split("\n").filter(Boolean)
    .map(l => JSON.parse(l));
}

// ── Call skills shim ──────────────────────────────────────────────────────────
// Contract: model + UNCHANGED user prompt + optional systemContext (here: a
// short metadata TEASER, not the full recipe — see buildSkillTeaser() below).
// skills_shim.js MUST route systemContext to the provider-specific
// system/developer instruction channel (Anthropic `system`, OpenAI
// `developer`/`system` message, Gemini/OpenRouter `system_instruction` — see
// skills_shim.js's file header). If the shim doesn't recognize the field it
// will likely just be dropped silently, which is exactly the failure mode
// the teaser smoke test in preflightChecks() is meant to catch.
//
// Payload is passed via a temp file, NOT embedded in the command string —
// same Windows ~8191-char command-line limit rationale as
// recipe_benchmark_runner.js. The teaser here is much shorter than a full
// recipe.md, but tool-result content (recipe.md's own text once fetched,
// plus any code/*.R files) accumulates in the conversation across turns, so
// the same safety margin applies.
function callBaseline(model, prompt, systemContext = null) {
  return new Promise((resolve) => {
    const body = systemContext != null ? { model, prompt, systemContext } : { model, prompt };
    const shimPath = path.join(__dirname, "scripts/skills_shim.js");

    fs.mkdirSync(CONFIG.results_dir, { recursive: true });
    const payloadPath = path.join(
      CONFIG.results_dir,
      `_payload_${Date.now()}_${Math.random().toString(36).slice(2, 8)}.json`
    );
    fs.writeFileSync(payloadPath, JSON.stringify(body), "utf8");

    const cleanup = () => { try { fs.unlinkSync(payloadPath); } catch {} };

    exec(
      `node "${shimPath}" --payload-file "${payloadPath}"`,
      // 300s, up from 180s — when systemContext is present, skills_shim.js
      // now runs a multi-turn tool-calling loop (up to 8 view() calls
      // round trips before a final answer), each turn a full API round trip
      // plus retry backoff on transient errors. 180s was already tight for a
      // single call under retry; it's not enough headroom for up to 6 calls.
      { timeout: 300000 },
      (error, stdout) => {
        cleanup();
        try {
          const parsed = JSON.parse(stdout?.trim().split("\n").filter(l => l.startsWith("{")).pop() || "{}");
          resolve(parsed);
        } catch {
          const timedOut = error?.killed || error?.signal;
          resolve({
            result: null,
            error: timedOut ? `SHIM_TIMEOUT: process killed after ${error?.signal || "timeout"}` : (error?.message || "parse error")
          });
        }
      }
    );
  });
}

// ── Recipe loading & injection ───────────────────────────────────────────────

function compactMarkdown(text) {
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
    chars: value.length,
    estimated_tokens: estimateTokens(value),
    lines: value.length === 0 ? 0 : value.split(/\r?\n/).length,
  };
}

// Loaded once per run — recipe.md does not change per-prompt or per-model.
// `raw` is the complete, unmodified file — used for two things: (1) hashing
// for provenance (recipeDoc.sha256, saved in results so you can prove which
// exact recipe revision a run used), and (2) extracting the short teaser
// below. Unlike recipe_benchmark_runner.js, `raw` itself is NEVER sent to a
// model directly here — models only ever see it by calling view("recipe.md")
// themselves, exactly like recipe_benchmark_runner.js's principle of "no
// summarisation", just applied one layer further back: we don't even
// pre-decide THAT the model reads it, only make it discoverable.
function loadRecipeDocument() {
  if (!fs.existsSync(CONFIG.recipe_file)) {
    throw new Error(`Recipe file not found: ${CONFIG.recipe_file}`);
  }
  const raw = fs.readFileSync(CONFIG.recipe_file, "utf8");
  const stats = textStats(raw);

  log(`  [recipe] recipe.md raw_chars=${raw.length} estimated_tokens=${stats.estimated_tokens} sha256=${sha256(raw).slice(0, 12)}… (NOTE: not injected directly in this script — see buildSkillTeaser)`);

  if (stats.chars > CONFIG.recipe_warning_chars) {
    log(`  [recipe] WARNING: recipe.md is ${stats.chars} chars, above warning threshold ${CONFIG.recipe_warning_chars}`);
  }
  if (stats.chars > CONFIG.recipe_max_chars) {
    throw new Error(`recipe.md is ${stats.chars} chars, above BENCH_RECIPE_MAX_CHARS=${CONFIG.recipe_max_chars}. This run refuses to silently truncate it — shorten recipe.md or raise the limit explicitly.`);
  }
  return { raw, stats, sha256: sha256(raw) };
}

// Pulls `id:` / `name:` out of recipe.md's YAML frontmatter (between the
// first two `---` lines) without a full YAML parser — the frontmatter here
// is flat key: value pairs, a regex is sufficient and avoids a dependency.
function extractFrontmatterField(raw, field) {
  const frontmatterMatch = raw.match(/^---\n([\s\S]*?)\n---/);
  if (!frontmatterMatch) return null;
  const line = frontmatterMatch[1].split("\n").find(l => l.trim().startsWith(`${field}:`));
  return line ? line.split(":").slice(1).join(":").trim() : null;
}

// Pulls the body paragraph under a given "## Heading" (up to the next "##").
function extractSection(raw, heading) {
  const re = new RegExp(`##\\s+${heading.replace(/[.*+?^${}()|[\\]\\\\]/g, "\\$&")}\\s*\\n([\\s\\S]*?)(?=\\n##\\s|$)`, "i");
  const match = raw.match(re);
  return match ? match[1].trim() : null;
}

// Builds the SHORT metadata-only teaser sent as systemContext — this is the
// entire point of this script versus recipe_benchmark_runner.js. It gives
// the model just enough to know a relevant skill exists and roughly what it
// covers (mirroring a real Skill's name + description, the only thing
// visible before selection), NOT the skill's actual instructions. The model
// must use the view tool to find and read the real thing itself. Note this
// deliberately does NOT reveal the filename "recipe.md" — the model has to
// list the folder and recognize the main instruction file itself, the same
// way discovering an unfamiliar skill folder actually works.
function buildSkillTeaser(recipeRaw) {
  const name = extractFrontmatterField(recipeRaw, "name") || "a domain-specific analysis skill";
  const whatItDoes = extractSection(recipeRaw, "What this recipe does") || "(description unavailable)";

  return [
    "A skill is available and has already been determined to be relevant to the task you're about to be given.",
    "",
    `Name: ${name}`,
    "",
    `What it does: ${whatItDoes}`,
    "",
    "You do not yet have this skill's full instructions. First use view(\".\") to see what's in the skill's folder, identify its main instruction file, and read it in full. Then read only the smallest number of supporting files needed for the SPECIFIC task you're given — usually one to three. Do not inspect every referenced pipeline script; only what this particular task actually requires. Once you have enough information, stop using tools and return the final executable R code.",
  ].join("\n");
}


// ── Task-specific validation (identical to benchmark_runner.js) ─────────────

function buildValidationR(validation = {}) {
  const checks = [];
  const explicitObj = validation.primary_object;
  const needsPrimaryObject =
    validation.min_cells != null ||
    validation.min_features != null ||
    validation.required_reductions?.length ||
    validation.require_clusters ||
    validation.required_columns?.length;

  if (needsPrimaryObject && !explicitObj) {
    checks.push(`
.bench_seurat_candidates <- Filter(function(x) {
  v <- tryCatch(get(x, envir = .GlobalEnv), error = function(e) NULL)
  inherits(v, "Seurat")
}, ls(envir = .GlobalEnv))
if (length(.bench_seurat_candidates) == 0) stop("VALIDATION_FAILED: no Seurat object found in environment")
.bench_obj <- get(.bench_seurat_candidates[[1]], envir = .GlobalEnv)
`.trim());
  }

  const obj = explicitObj || ".bench_obj";

  if (validation.required_objects?.length) {
    for (const o of validation.required_objects) {
      checks.push(`if (!exists("${o}")) stop("VALIDATION_FAILED: required object '${o}' not found")`);
    }
  }
  if (validation.required_classes) {
    for (const [o, cls] of Object.entries(validation.required_classes)) {
      checks.push(`if (!inherits(${o}, "${cls}")) stop("VALIDATION_FAILED: '${o}' is not class ${cls}")`);
    }
  }
  if (validation.min_cells != null) {
    checks.push(`if (ncol(${obj}) < ${validation.min_cells}) stop("VALIDATION_FAILED: too few cells")`);
  }
  if (validation.min_features != null) {
    checks.push(`if (nrow(${obj}) < ${validation.min_features}) stop("VALIDATION_FAILED: too few features")`);
  }
  if (validation.required_reductions?.length) {
    for (const r of validation.required_reductions) {
      checks.push(`if (!"${r}" %in% names(${obj}@reductions)) stop("VALIDATION_FAILED: missing reduction '${r}'")`);
    }
  }
  if (validation.require_clusters) {
    const minClusters = validation.min_clusters || 2;
    checks.push(`if (length(unique(Idents(${obj}))) < ${minClusters}) stop("VALIDATION_FAILED: insufficient clusters")`);
  }
  if (validation.required_columns?.length) {
    for (const col of validation.required_columns) {
      checks.push(`if (!"${col}" %in% colnames(${obj}@meta.data)) stop("VALIDATION_FAILED: missing metadata column '${col}'")`);
    }
  }

  return checks.join("\n");
}

// ── Execute R code (sandboxed — identical isolation model to benchmark_runner.js) ──

const PACKAGE_CHECK_PREFIX = `
required_packages <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr","readr","ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(paste0("BENCH_ENV_ERROR: missing required packages: ", paste(missing_packages, collapse = ", "),
              ". Run setup_env.R inside the micromamba environment first."))
}
`.trim();

function executeRCode(rCode, datasetDir, validation = {}) {
  const tmpId   = Math.random().toString(36).slice(2, 8);
  const tmpDir  = path.join(CONFIG.results_dir, `rexec_recipe_${tmpId}`);
  fs.mkdirSync(tmpDir, { recursive: true });
  const scriptPath = path.join(tmpDir, "bench.R");

  rCode = rCode
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']\.\/hg19\/?["']\s*\)/gi, 'Read10X(data.dir = ".")')
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']hg19\/?["']\s*\)/gi,    'Read10X(data.dir = ".")')
    .replace(/(data\.dir\s*=\s*)["'][^"']*(?:path|your|data)[^"']*["']/gi, '$1"."')
    .replace(/(data_dir\s*<-\s*)["'][^"']*(?:path|your|data)[^"']+["']/gi, '$1"."');

  const validationCode = buildValidationR(validation);
  const hasValidation   = validationCode.length > 0;

  const fullScript = [
    `setwd("${datasetDir.replace(/\\/g, "/")}")`,
    PACKAGE_CHECK_PREFIX,
    `tryCatch({`,
    `  suppressMessages(suppressWarnings({`,
    rCode,
    `  }))`,
    hasValidation ? validationCode : "",
    hasValidation ? `  cat("BENCH_VALIDATION_OK\\n")` : "",
    `  cat("BENCH_EXEC_OK\\n")`,
    `}, error = function(e) {`,
    `  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\\n", file=stderr())`,
    `})`,
  ].join("\n");

  fs.writeFileSync(scriptPath, fullScript);
  log(`  [R script] ${scriptPath}`);
  log(`  [R exec] "${CONFIG.r_exec}" "${scriptPath}"`);

  return new Promise((resolve) => {
    const start = Date.now();
    execFile(CONFIG.r_exec, ["--vanilla", scriptPath], { timeout: CONFIG.timeout_ms },
      (error, stdout, stderr) => {
        const duration_ms      = Date.now() - start;
        const executionSuccess = !error && (stdout || "").includes("BENCH_EXEC_OK");
        const validationSuccess = hasValidation ? (stdout || "").includes("BENCH_VALIDATION_OK") : true;
        const success = executionSuccess && validationSuccess;
        const failureType = success ? null : classifyRFailure(stderr, stdout, error);

        if (!success) {
          log(`  [R exit] error=${error?.code} duration=${duration_ms}ms failureType=${failureType}`);
          log(`  [R stderr] ${tail(stderr, 600)}`);
          log(`  [R stdout] ${tail(stdout, 400)}`);
        } else {
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        resolve({
          success,
          executionSuccess,
          validationSuccess,
          validated: hasValidation,
          failureType,
          stderr: stderr || "",
          duration_ms,
        });
      }
    );
  });
}

function classifyRFailure(stderr = "", stdout = "", error = null) {
  const combined = `${stderr}\n${stdout}`;
  if (/BENCH_ENV_ERROR:/i.test(combined))                                  return "environment_error";
  if (/VALIDATION_FAILED:/i.test(combined))                                return "validation_mismatch";
  if (error?.killed || error?.signal)                                      return "timeout_or_termination";
  if (/duplicate rownames|duplicate 'row\.names'/i.test(combined))         return "duplicate_feature_names";
  if (/package .* is required|there is no package called/i.test(combined))return "missing_package";
  if (/could not find function/i.test(combined))                          return "missing_function";
  if (/object ['"]?[\w.]+['"]? not found/i.test(combined))                 return "missing_object";
  if (/directory provided does not exist|cannot open the connection/i.test(combined)) return "input_path_error";
  if (!stdout && !stderr && error)                                        return "no_output_captured";
  return "r_execution_error";
}

function tail(text, n) {
  if (!text) return "(empty)";
  return text.length > n ? "…" + text.slice(-n) : text;
}

// ── Metrics (identical to benchmark_runner.js — same grading = valid diff) ──

const R_REGISTRY_PACKAGES = [
  "base","stats","utils","graphics","grDevices","methods",
  "Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr",
  "readr","ggrepel","DoubletFinder","scDblFinder","SingleCellExperiment"
];

function getRFunctionRegistry() {
  const registryScript = `
packages <- c(${R_REGISTRY_PACKAGES.map(p => `"${p}"`).join(", ")})
all_functions <- unique(unlist(lapply(packages, function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) { message("REGISTRY_SKIP: ", pkg, " not installed"); return(character()) }
  ns <- asNamespace(pkg)
  exports <- getNamespaceExports(pkg)
  exports[vapply(exports, function(fn) {
    obj <- tryCatch(get(fn, envir = ns), error = function(e) NULL)
    is.function(obj)
  }, logical(1))]
})))
cat("REGISTRY_COUNT:", length(all_functions), "\\n", file = stderr())
cat(paste(all_functions, collapse = "\\n"))
`.trim();

  const tmpPath = path.join(CONFIG.results_dir, `_registry_query_recipe_${Date.now()}.R`);
  fs.mkdirSync(CONFIG.results_dir, { recursive: true });
  fs.writeFileSync(tmpPath, registryScript);

  try {
    const output = execFileSync(CONFIG.r_exec, ["--vanilla", tmpPath], { encoding: "utf8", timeout: 120000 });
    const entries = output.split("\n").map(l => l.trim()).filter(Boolean);
    const set = new Set(entries);
    log(`  R function registry built: ${set.size} exported functions from ${R_REGISTRY_PACKAGES.length} packages`);
    if (set.size < 100) {
      log(`  ⚠️  Registry looks suspiciously small (${set.size})`);
    }
    return set;
  } catch (e) {
    log(`  ⚠️  Could not build R function registry: ${e.message}`);
    return null;
  } finally {
    try { fs.unlinkSync(tmpPath); } catch {}
  }
}

const NON_FUNCTION_KEYWORDS = new Set([
  "if","else","for","while","repeat","function","return","switch","next","break"
]);

function extractDefinedFunctions(code) {
  const defined = new Set();
  const pattern = /\b([A-Za-z][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(/g;
  for (const m of code.matchAll(pattern)) defined.add(m[1]);
  return defined;
}

function stripRComments(code) {
  return code.split("\n").map(line => line.replace(/#.*$/, "")).join("\n");
}

function calcGroundingRate(code, registry) {
  if (!code) return { rate: 0, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: !!registry };
  if (!registry) return { rate: null, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: false };

  code = stripRComments(code);
  const definedFunctions = extractDefinedFunctions(code);
  const calls = [...code.matchAll(/\b([A-Za-z][A-Za-z0-9_.]*)\s*\(/g)]
    .map(m => m[1])
    .filter(name => !NON_FUNCTION_KEYWORDS.has(name));

  if (!calls.length) return { rate: 1, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: true };

  const grounded = [], ungrounded = [];
  for (const name of calls) {
    if (registry.has(name) || definedFunctions.has(name)) grounded.push(name);
    else ungrounded.push(name);
  }

  return {
    rate: grounded.length / calls.length,
    groundedCalls: grounded,
    ungroundedCalls: [...new Set(ungrounded)],
    totalCalls: calls.length,
    registryAvailable: true,
  };
}

const PIPELINE_STEPS = [
  { name: "data_loading",     patterns: ["Read10X","ReadMtx","ReadMTX","readMM","read.table","read.delim","read.csv","read_tsv","readLines"] },
  { name: "seurat_object",    patterns: ["CreateSeuratObject"] },
  { name: "qc",               patterns: ["PercentageFeatureSet","VlnPlot","FeatureScatter","subset"] },
  { name: "normalisation",    patterns: ["NormalizeData","SCTransform"] },
  { name: "variable_feats",   patterns: ["FindVariableFeatures"] },
  { name: "scaling",          patterns: ["ScaleData"] },
  { name: "pca",              patterns: ["RunPCA"] },
  { name: "graph",            patterns: ["FindNeighbors"] },
  { name: "clustering",       patterns: ["FindClusters"] },
  { name: "embedding",        patterns: ["RunUMAP","RunTSNE"] },
  { name: "marker_detection", patterns: ["FindAllMarkers","FindMarkers"] },
];
const STEP_MAP = new Map(PIPELINE_STEPS.map(s => [s.name, s]));

function findFirstPatternPosition(code, patterns) {
  const positions = patterns.map(p => code.indexOf(p)).filter(p => p >= 0);
  return positions.length ? Math.min(...positions) : -1;
}

function calcWorkflowScore(code, expectedWorkflow = []) {
  if (!expectedWorkflow.length) return { score: null, presenceScore: null, orderScore: null, steps: [] };
  if (!code) return { score: 0, presenceScore: 0, orderScore: 0, steps: [] };

  code = stripRComments(code);

  const evaluated = expectedWorkflow.map(stepName => {
    const step = STEP_MAP.get(stepName);
    if (!step) return { name: stepName, present: false, position: -1, unknown: true };
    const position = findFirstPatternPosition(code, step.patterns);
    return { name: stepName, present: position >= 0, position };
  });

  const presenceScore = evaluated.filter(s => s.present).length / evaluated.length;
  const presentSteps  = evaluated.filter(s => s.present);

  let correctPairs = 0, totalPairs = 0;
  for (let i = 0; i < presentSteps.length - 1; i++) {
    totalPairs++;
    if (presentSteps[i].position < presentSteps[i + 1].position) correctPairs++;
  }
  const orderScore = totalPairs === 0 ? presenceScore : correctPairs / totalPairs;

  return { score: 0.7 * presenceScore + 0.3 * orderScore, presenceScore, orderScore, steps: evaluated };
}

function passAtK(n, c, k) {
  if (n - c < k) return 1.0;
  let num = 1, den = 1;
  for (let i = 0; i < k; i++) {
    num *= (n - c - i);
    den *= (n - i);
  }
  return 1 - num / den;
}

function normalizeCodeForComparison(code) {
  return (code || "").replace(/#.*$/gm, "").replace(/\s+/g, " ").trim();
}

function countUniqueGenerations(attempts) {
  const unique = new Set(attempts.map(a => normalizeCodeForComparison(a.code)).filter(Boolean));
  return unique.size;
}

function extractRFromText(text) {
  if (!text) return null;
  text = text.replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "").trim();
  let match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  match = text.match(/```[a-z]*\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) return text;
  return null;
}

// ── Pass@k suite (identical grading logic to benchmark_runner.js; the only
//    difference is upstream — prompts arrive with the recipe attached as
//    p.system_context, while p.prompt is left byte-identical to the vanilla
//    baseline. See buildSkillTeaser() and the main flow below.) ─────

async function runPassK(modelName, callerFn, prompts, k) {
  log(`\n[Pass@k] model=${modelName}  k=${k}  prompts=${prompts.length}  condition=skills_progressive_disclosure (teaser + view tool)`);
  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < k; i++) {
      // p.prompt is the ORIGINAL, unmodified task text (identical to the
      // vanilla-baseline run). p.system_context carries the recipe and is
      // routed separately by callerFn/callBaseline/the shim.
      const res  = await callerFn(p.prompt, p.system_context);
      const code = res.alreadyExtracted
        ? (res.result || null)
        : extractRFromText(res.result || "");
      const toolCallsMade = res.tool_calls_made || [];
      const toolBudgetExhausted = !!res.tool_budget_exhausted;
      if (toolCallsMade.length > 0) {
        log(`  [${p.id}] tool calls: ${toolCallsMade.map(c => `${c?.path}(${c?.result})`).join(", ")}${toolBudgetExhausted ? " [forced final synthesis]" : ""}`);
      }
      if (!code) {
        // Distinct failure_type when the exploration budget was exhausted
        // AND the model still couldn't produce extractable code even on the
        // forced, tools-disabled final turn — a genuinely different (and
        // with the forced-synthesis fix, now much rarer) failure mode than
        // ordinary extraction failure. Keeps "ran out of budget but still
        // couldn't answer" distinguishable from "produced malformed output
        // for unrelated reasons" in later analysis.
        const failureType = toolBudgetExhausted ? "tool_budget_exhausted_no_code" : "no_code_extracted";
        log(`  ⚠️  [${p.id}] ${failureType} — model_error=${res.error || "(none)"} raw_preview=${JSON.stringify((res.raw_preview || "(not captured)").slice(0, 300))}`);
        attempts.push({ passed: false, error: "no_r_code_extracted", failure_type: failureType, model_error: res.error || null, raw_preview: res.raw_preview || null, tool_calls_made: toolCallsMade, tool_budget_exhausted: toolBudgetExhausted });
        continue;
      }
      const exec = await executeRCode(code, CONFIG.dataset_dir, p.validation || {});
      attempts.push({
        passed:            exec.success,
        code,
        executionSuccess:  exec.executionSuccess,
        validationSuccess: exec.validationSuccess,
        validated:         exec.validated,
        failure_type:      exec.failureType,
        error:             exec.success ? null : (exec.stderr || "unknown error").trim(),
        duration_ms:       exec.duration_ms,
        tool_calls_made:   toolCallsMade,
        tool_budget_exhausted: toolBudgetExhausted,
      });
    }
    const c = attempts.filter(a => a.passed).length;
    const score = passAtK(k, c, k);
    const uniqueGenerations = countUniqueGenerations(attempts);
    log(`  [${p.id}] pass@${k}=${score.toFixed(3)}  (${c}/${k} passed, ${uniqueGenerations}/${k} unique generations)`);
    results.push({
      id: p.id,
      score,
      attempt_success_rate: attempts.length ? c / attempts.length : 0,
      unique_generations: uniqueGenerations,
      expected_workflow: p.expected_workflow || [],
      validated: attempts.some(a => a.validated),
      attempts,
    });
  }

  const avg = results.reduce((s, r) => s + r.score, 0) / results.length;
  const avgAttemptSuccessRate = results.reduce((s, r) => s + r.attempt_success_rate, 0) / results.length;
  log(`  avg Pass@${k} = ${avg.toFixed(3)}   avg attempt success rate = ${avgAttemptSuccessRate.toFixed(3)}`);
  return { avg, avgAttemptSuccessRate, results };
}

function runCodeQuality(modelName, passKResults, registry) {
  log(`\n[CodeQuality] model=${modelName} (teaser + view-tool attempts, registry ${registry ? "available" : "UNAVAILABLE"})`);
  const scores = [];

  for (const r of passKResults) {
    const attemptScores = r.attempts.map((attempt, i) => {
      const code = attempt.code || null;
      const workflowResult  = code ? calcWorkflowScore(code, r.expected_workflow) : { score: 0, presenceScore: 0, orderScore: 0, steps: [] };
      const groundingResult = code ? calcGroundingRate(code, registry)            : { rate: 0, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: !!registry };
      return { attempt: i + 1, passed: attempt.passed, workflow: workflowResult, grounding: groundingResult };
    });

    const workflowVals  = attemptScores.map(a => a.workflow.score).filter(v => v !== null);
    const groundingVals = attemptScores.map(a => a.grounding.rate).filter(v => v !== null);
    const avgWorkflow    = workflowVals.length  ? workflowVals.reduce((s, v) => s + v, 0) / workflowVals.length : null;
    const avgGrounding   = groundingVals.length ? groundingVals.reduce((s, v) => s + v, 0) / groundingVals.length : null;

    const successfulOnly = attemptScores.filter(a => a.passed);
    const successfulWorkflowVals  = successfulOnly.map(a => a.workflow.score).filter(v => v !== null);
    const successfulGroundingVals = successfulOnly.map(a => a.grounding.rate).filter(v => v !== null);
    const successfulWorkflow  = successfulWorkflowVals.length  ? successfulWorkflowVals.reduce((s, v) => s + v, 0) / successfulWorkflowVals.length : null;
    const successfulGrounding = successfulGroundingVals.length ? successfulGroundingVals.reduce((s, v) => s + v, 0) / successfulGroundingVals.length : null;

    scores.push({ id: r.id, workflow: avgWorkflow, grounding: avgGrounding, successful_workflow: successfulWorkflow, successful_grounding: successfulGrounding, attempts: attemptScores });
  }

  const validWorkflow  = scores.map(s => s.workflow).filter(v => v !== null);
  const validGrounding = scores.map(s => s.grounding).filter(v => v !== null);
  const avgWorkflow  = validWorkflow.length  ? validWorkflow.reduce((s, v) => s + v, 0) / validWorkflow.length : null;
  const avgGrounding = validGrounding.length ? validGrounding.reduce((s, v) => s + v, 0) / validGrounding.length : null;

  log(`  avg workflow success rate: ${avgWorkflow !== null ? avgWorkflow.toFixed(3) : "N/A"}`);
  log(`  avg grounding rate:        ${avgGrounding !== null ? avgGrounding.toFixed(3) : "N/A"}`);
  return { avgWorkflow, avgGrounding, scores };
}

// ── Model registry ────────────────────────────────────────────────────────────
// Deliberately EXCLUDES rolv / rolv-nopreexec: those already perform their own
// RAG-based context retrieval, so prepending the full recipe on top would
// conflate two different context-injection mechanisms in one condition.
//
// Also EXCLUDES rolv-raw: excluded from recipe_benchmark_runner.js after its
// dry run consistently failed with "fetch failed" (ROLV endpoint
// unreachable at test time), and never re-added here either — skills_shim.js
// has no callRolvRaw() implementation at all (unlike baseline_shim.js, which
// kept a working-but-unused one). Add ROLV raw support to skills_shim.js
// from scratch if/when ROLV connectivity and view-tool-via-ROLV support are
// separately confirmed.

const API_MODELS = {
  groq:                (p, sc) => callBaseline("groq", p, sc),
  "gpt-5.5":           (p, sc) => callBaseline("gpt-5.5", p, sc),
  "gpt-5.4-mini":      (p, sc) => callBaseline("gpt-5.4-mini", p, sc),
  "gemini-openrouter": (p, sc) => callBaseline("gemini-openrouter", p, sc),
  "gemini-flash-lite": (p, sc) => callBaseline("gemini-flash-lite", p, sc),
  "claude-opus-4-8":   (p, sc) => callBaseline("claude-opus-4-8", p, sc),
  "claude-sonnet-5":   (p, sc) => callBaseline("claude-sonnet-5", p, sc),
  "claude-haiku-4-5":  (p, sc) => callBaseline("claude-haiku-4-5", p, sc),
};

const MODEL_DISPLAY_NAMES = {
  "groq":              "GPT-OSS 20B via Groq + skill (view tool)",
  "gpt-5.5":           "GPT-5.5 (OpenAI flagship) + skill (view tool)",
  "gpt-5.4-mini":      "GPT-5.4 Mini (OpenAI) + skill (view tool)",
  "gemini-openrouter": "Gemini 3.5 Flash (via OpenRouter) + skill (view tool)",
  "gemini-flash-lite": "Gemini 3.1 Flash-Lite (via OpenRouter) + skill (view tool)",
  "claude-opus-4-8":   "Claude Opus 4.8 (Anthropic direct) + skill (view tool)",
  "claude-sonnet-5":   "Claude Sonnet 5 (Anthropic direct) + skill (view tool)",
  "claude-haiku-4-5":  "Claude Haiku 4.5 (Anthropic direct) + skill (view tool)",
};

const MODEL_REQUIRED_KEYS = {
  "groq":              ["GROQ_API_KEY"],
  "gpt-5.5":           ["OPENAI_API_KEY"],
  "gpt-5.4-mini":      ["OPENAI_API_KEY"],
  "gemini-openrouter": ["OPENROUTER_API_KEY"],
  "gemini-flash-lite": ["OPENROUTER_API_KEY"],
  "claude-opus-4-8":   ["ANTHROPIC_API_KEY"],
  "claude-sonnet-5":   ["ANTHROPIC_API_KEY"],
  "claude-haiku-4-5":  ["ANTHROPIC_API_KEY"],
};

const KEY_LABELS = {
  ROLV_API_URL:       "ROLV API URL",
  ROLV_JWT_TOKEN:     "ROLV JWT token",
  OPENAI_API_KEY:     "OpenAI key",
  OPENROUTER_API_KEY: "OpenRouter key",
  GROQ_API_KEY:       "Groq key",
  ANTHROPIC_API_KEY:  "Anthropic key",
};

// ── Preflight checks ──────────────────────────────────────────────────────────

async function preflightChecks() {
  log("═".repeat(60));
  log("PREFLIGHT CHECKS (skills progressive-disclosure run)");
  log("═".repeat(60));

  const errors = [], warnings = [];

  const selectedModels = CONFIG.model === "all"
    ? Object.keys(MODEL_REQUIRED_KEYS)
    : CONFIG.model.split(",").map(m => m.trim());

  log(`  Models in scope: ${selectedModels.join(", ")}`);

  const unknownModels = selectedModels.filter(m => !Object.prototype.hasOwnProperty.call(MODEL_REQUIRED_KEYS, m));
  if (unknownModels.length) {
    errors.push(`Unknown model key(s): ${unknownModels.join(", ")} — check spelling against: ${Object.keys(MODEL_REQUIRED_KEYS).join(", ")}`);
  }

  const neededKeys = new Set();
  for (const m of selectedModels) for (const key of (MODEL_REQUIRED_KEYS[m] || [])) neededKeys.add(key);
  for (const key of neededKeys) {
    const val = process.env[key];
    if (!val || val.startsWith("sk-replace") || val.startsWith("sk-or-replace")) {
      errors.push(`Missing or placeholder: ${key} (${KEY_LABELS[key] || key})`);
    } else {
      log(`  ✅ ${key} present`);
    }
  }

  // Recipe file
  try {
    loadRecipeDocument();
    log(`  ✅ recipe.md loaded from ${CONFIG.recipe_file}`);
  } catch (e) {
    errors.push(e.message);
  }

  // Baseline file for comparison (warning only — run can still proceed)
  if (!CONFIG.skipComparison) {
    if (!fs.existsSync(CONFIG.baseline_file)) {
      warnings.push(`Baseline file not found: ${CONFIG.baseline_file} — comparison report will be skipped (pass --baseline <path> or --skip-comparison to silence this)`);
    } else {
      log(`  ✅ Baseline file found: ${CONFIG.baseline_file}`);
    }
  }

  // R reachable + packages
  let rVersionOk = false;
  try {
    const versionOut = execFileSync(CONFIG.r_exec, ["--version"], { encoding: "utf8", timeout: 15000 });
    log(`  ✅ R reachable: ${versionOut.split("\n")[0].trim()}`);
    rVersionOk = true;
  } catch (e) {
    errors.push(`R binary not reachable via "${CONFIG.r_exec}": ${e.message}`);
  }

  if (rVersionOk) {
    try {
      const pkgCheckScript = `
required_packages <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr","readr","ggrepel","DoubletFinder","scDblFinder","SingleCellExperiment")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) cat(pkg, ":", as.character(packageVersion(pkg)), "\\n") else cat(pkg, ": MISSING\\n")
}
`.trim();
      const pkgTmpPath = path.join(CONFIG.results_dir, `_pkgcheck_recipe_${Date.now()}.R`);
      fs.mkdirSync(CONFIG.results_dir, { recursive: true });
      fs.writeFileSync(pkgTmpPath, pkgCheckScript);
      let pkgOut;
      try {
        // 60s, up from 30s — this spawns a fresh Rscript process that loads
        // Seurat + 12 other packages from scratch, which is legitimately
        // slow (several seconds) and gets slower under system load (e.g.
        // right after several back-to-back benchmark runs each spawning
        // their own R processes). A timeout here does NOT mean R itself is
        // unreachable — rVersionOk is already true by this point — it means
        // this specific package-loading check ran out of time, which is
        // usually transient load, not a real environment problem.
        pkgOut = execFileSync(CONFIG.r_exec, ["--vanilla", pkgTmpPath], { encoding: "utf8", timeout: 60000 });
      } finally {
        try { fs.unlinkSync(pkgTmpPath); } catch {}
      }
      const missing = [];
      for (const line of pkgOut.split("\n").map(l => l.trim()).filter(Boolean)) {
        if (line.endsWith("MISSING")) missing.push(line.split(":")[0].trim());
      }
      if (missing.length) errors.push(`Missing R packages: ${missing.join(", ")}`);
    } catch (e) {
      const timedOut = e?.signal || /ETIMEDOUT/i.test(e?.message || "");
      errors.push(
        timedOut
          ? `R package-availability check timed out (R itself IS reachable — see the ✅ line above). This is usually transient system load, e.g. right after several back-to-back R process spawns, not a real problem with your R/Seurat install. Try again in a few seconds, ideally not immediately after another full run: ${e.message}`
          : `R package-availability check failed: ${e.message}`
      );
    }
  }

  // Dataset dir
  if (!fs.existsSync(CONFIG.dataset_dir)) {
    errors.push(`Dataset directory not found: ${CONFIG.dataset_dir}`);
  } else {
    log(`  ✅ Dataset directory found: ${CONFIG.dataset_dir}`);
  }

  // Prompts file + validation coverage
  if (!fs.existsSync(CONFIG.prompts_file)) {
    errors.push(`Prompts file not found: ${CONFIG.prompts_file}`);
  } else {
    log(`  ✅ Prompts file found`);
    try {
      const promptsPreview = loadJSONL(CONFIG.prompts_file);
      const missingValidation = promptsPreview.filter(p => !p.validation || Object.keys(p.validation).length === 0);
      const missingWorkflow   = promptsPreview.filter(p => !p.expected_workflow || p.expected_workflow.length === 0);
      const bucket = CONFIG.allowPartialPrompts ? warnings : errors;
      if (missingValidation.length) {
        bucket.push(`${missingValidation.length}/${promptsPreview.length} prompts have no 'validation' field: ${missingValidation.map(p => p.id).join(", ")}` + (CONFIG.allowPartialPrompts ? "" : " (pass --allow-partial-prompts to run anyway)"));
      }
      if (missingWorkflow.length) {
        bucket.push(`${missingWorkflow.length}/${promptsPreview.length} prompts have no 'expected_workflow': ${missingWorkflow.map(p => p.id).join(", ")}` + (CONFIG.allowPartialPrompts ? "" : " (pass --allow-partial-prompts to run anyway)"));
      }
      const knownSteps = new Set(PIPELINE_STEPS.map(s => s.name));
      const badSteps = [];
      for (const p of promptsPreview) for (const step of (p.expected_workflow || [])) if (!knownSteps.has(step)) badSteps.push(`${p.id}:"${step}"`);
      if (badSteps.length) errors.push(`Unknown expected_workflow step name(s): ${badSteps.join(", ")}`);
    } catch (e) {
      warnings.push(`Could not pre-check prompt validation coverage: ${e.message}`);
    }
  }

  // Live smoke tests
  if (!errors.length) {
    log("\n  Running live API smoke tests (1 call each, selected models only)...");
    // Three separate mechanisms are being verified, and they can fail
    // independently: (1) systemContext (the teaser) actually reaching the
    // model's system/developer channel, (2) the view tool actually being
    // usable by that provider/route for BOTH listing a directory and
    // reading a file — not just one or the other, since a provider could
    // technically call the tool with a directory path and never try a file
    // path, or vice versa, and (3) content fetched via the tool actually
    // making it into the final answer, not just being requested and ignored.
    //
    // The prompt explicitly instructs both calls rather than hoping the
    // model spontaneously explores — the smoke test's synthetic teaser
    // (below) is deliberately tiny and doesn't describe the real skill, so
    // nothing would organically prompt exploration without being told to.
    // Real prompts DO get the real teaser (which explicitly tells the model
    // to explore), so real usage is expected to happen more naturally there
    // — this smoke test only proves the MECHANISM works end-to-end.
    //
    // Reference tag framing (not "secret") avoids the refusal behavior
    // GPT-5.5 showed for "secret code word" framing in an earlier version
    // of recipe_benchmark_runner.js's equivalent test — that was actually
    // proof systemContext WAS delivered and read, the model just declined
    // to leak content framed as confidential.
    const SMOKE_SYSTEM = "Background reference tag for this session: PINEAPPLE-42.";
    const SMOKE_PROMPT = 'First, call view(".") to list the current folder. Then call view("code/05_doublets.R") to read that file and find the value of pK_fixed defined there. Then write a single line of R code that assigns a named list to a variable called x containing two things: the reference tag mentioned in your instructions, and the pK_fixed value you just looked up. Example: x <- list(tag = "TAG", pk = 0.0)';
    for (const name of selectedModels) {
      if (!API_MODELS[name]) continue;
      try {
        const res = await API_MODELS[name](SMOKE_PROMPT, SMOKE_SYSTEM);
        const hasResult = res.result != null && res.result !== "";
        if (res.error && !hasResult) throw new Error(res.error);

        const calls = res.tool_calls_made || [];
        const echoedTag  = hasResult && /PINEAPPLE-?42/i.test(res.result);
        const echoedPk   = hasResult && /0\.01/.test(res.result);
        const usedList   = calls.some(c => c?.path === "." || c?.path === "./");
        const usedRead   = calls.some(c => c?.path === "code/05_doublets.R");

        if (hasResult && !echoedTag) {
          errors.push(`${name}: systemContext (teaser) was NOT confirmed applied (extracted R code did not contain the smoke-test reference tag) — skills_shim.js is likely not routing systemContext to this provider's system/developer field. Response preview: ${JSON.stringify((res.result || "").slice(0, 150))}`);
          log(`  ❌ ${name}: systemContext NOT confirmed applied`);
        } else if (!usedList || !usedRead) {
          const missing = [!usedList && "list", !usedRead && "read"].filter(Boolean).join(" and ");
          errors.push(`${name}: view tool did NOT complete a ${missing} call even though explicitly instructed to (tool_calls_made=${JSON.stringify(calls)}) — this provider/route likely doesn't support function-calling the way skills_shim.js sends it, or only handles one of directory-listing/file-reading. Response preview: ${JSON.stringify((res.result || "").slice(0, 150))}`);
          log(`  ❌ ${name}: view tool ${missing} call NOT confirmed`);
        } else if (!echoedPk) {
          errors.push(`${name}: view WAS called correctly (list + read) but the returned pK_fixed value (0.01) never appeared in the final answer — tool result may not be getting appended to the conversation correctly. Response preview: ${JSON.stringify((res.result || "").slice(0, 150))}`);
          log(`  ❌ ${name}: tool called but result not reflected in final answer`);
        } else {
          log(`  ✅ ${name} reachable, teaser confirmed applied, view tool confirmed working (list + read)`);
        }
      } catch (e) {
        errors.push(`${name} smoke test FAILED: ${e.message}`);
        log(`  ❌ ${name}: ${e.message}`);
      }
    }
  }

  log("\n" + "─".repeat(60));
  if (errors.length) {
    log("❌ PREFLIGHT FAILED — fix these before running:");
    errors.forEach(e => log(`   • ${e}`));
    log("\nAborting. No prompts were run.");
    process.exit(1);
  }
  if (warnings.length) warnings.forEach(w => log(`  ⚠️  ${w}`));
  log("✅ All preflight checks passed — starting skills progressive-disclosure benchmark run.");
  log("─".repeat(60));
}

// ── Comparison against frozen baseline ───────────────────────────────────────

function pct(v) { return v === null || v === undefined ? "N/A" : (v * 100).toFixed(1) + "%"; }

function buildComparison(recipeResults) {
  if (CONFIG.skipComparison) return null;
  if (!fs.existsSync(CONFIG.baseline_file)) {
    log(`\n⚠️  Skipping comparison — baseline file not found: ${CONFIG.baseline_file}`);
    return null;
  }

  let baseline;
  try {
    baseline = JSON.parse(fs.readFileSync(CONFIG.baseline_file, "utf8"));
  } catch (e) {
    log(`\n⚠️  Skipping comparison — could not parse baseline file: ${e.message}`);
    return null;
  }

  // Pass@k is NOT comparable across different k. Pass@3 is mathematically
  // >= Pass@1 for the same underlying per-attempt success rate (more
  // attempts raise the chance of at least one success) — so a k mismatch
  // alone would make the recipe condition look worse than it is, for a
  // reason unrelated to the recipe. Never fail silently on this: annotate
  // the comparison object AND warn loudly so it can't be misread later,
  // rather than blocking the run outright (k=1 dry-run-style comparisons
  // against a k=3 baseline are still useful for spotting gross breakage,
  // just not for a real Pass@k delta).
  const baselineK = baseline.config?.k;
  const kMismatch = baselineK != null && baselineK !== CONFIG.k;
  if (kMismatch) {
    log(`\n⚠️  K MISMATCH: this run used k=${CONFIG.k}, baseline used k=${baselineK}. Pass@k deltas below are NOT a like-for-like comparison — Pass@${baselineK} is structurally >= Pass@${CONFIG.k} for the same success rate, independent of any recipe effect. Treat deltas as directional/sanity-check only until re-run at matching k.`);
  }

  const rows = [];
  for (const [modelName, recipeRes] of Object.entries(recipeResults)) {
    const baseRes = baseline.results?.[modelName];
    if (!baseRes?.passK) {
      log(`  [compare] no baseline entry for model="${modelName}" — skipping`);
      continue;
    }
    const basePassK   = baseRes.passK.avg;
    const baseASR      = baseRes.passK.avgAttemptSuccessRate;
    const baseWorkflow = baseRes.quality?.avgWorkflow ?? null;
    const baseGrounding= baseRes.quality?.avgGrounding ?? null;

    const recipePassK   = recipeRes.passK.avg;
    const recipeASR      = recipeRes.passK.avgAttemptSuccessRate;
    const recipeWorkflow = recipeRes.quality?.avgWorkflow ?? null;
    const recipeGrounding= recipeRes.quality?.avgGrounding ?? null;

    rows.push({
      model: modelName,
      baseline_pass_k: basePassK,
      recipe_pass_k: recipePassK,
      delta_pass_k: recipePassK - basePassK,
      baseline_attempt_success_rate: baseASR,
      recipe_attempt_success_rate: recipeASR,
      delta_attempt_success_rate: recipeASR - baseASR,
      baseline_workflow: baseWorkflow,
      recipe_workflow: recipeWorkflow,
      delta_workflow: (baseWorkflow != null && recipeWorkflow != null) ? recipeWorkflow - baseWorkflow : null,
      baseline_grounding: baseGrounding,
      recipe_grounding: recipeGrounding,
      delta_grounding: (baseGrounding != null && recipeGrounding != null) ? recipeGrounding - baseGrounding : null,
      // Per-prompt breakdown, for spotting which specific tasks the recipe helped/hurt.
      per_prompt: (recipeRes.passK.results || []).map(rp => {
        const bp = (baseRes.passK.results || []).find(x => x.id === rp.id);
        return {
          id: rp.id,
          baseline_score: bp ? bp.score : null,
          recipe_score: rp.score,
          delta: bp ? rp.score - bp.score : null,
        };
      }),
    });
  }

  return {
    baseline_run_id: baseline.run_id,
    baseline_file: CONFIG.baseline_file,
    recipe_k: CONFIG.k,
    baseline_k: baselineK ?? null,
    k_mismatch: kMismatch,
    rows,
  };
}

function logComparisonTable(comparison) {
  if (!comparison) return;
  log("\n" + "─".repeat(68));
  log(`COMPARISON vs BASELINE (run_id=${comparison.baseline_run_id})`);
  if (comparison.k_mismatch) {
    log(`⚠️  K MISMATCH: this run k=${comparison.recipe_k} vs baseline k=${comparison.baseline_k} — deltas below are NOT a valid Pass@k comparison, directional only.`);
  }
  log("─".repeat(68));
  log(`${"Model".padEnd(30)} Baseline    +Recipe     Δ Pass@k`);
  log("─".repeat(68));
  for (const r of comparison.rows) {
    const label = (MODEL_DISPLAY_NAMES[r.model] || r.model).padEnd(30);
    const deltaStr = (r.delta_pass_k >= 0 ? "+" : "") + (r.delta_pass_k * 100).toFixed(1) + "pp";
    log(`${label} ${pct(r.baseline_pass_k).padEnd(11)} ${pct(r.recipe_pass_k).padEnd(11)} ${deltaStr}`);
  }
  log("─".repeat(68));
  log("(pp = percentage points. Positive Δ = recipe injection helped this model on this metric.)");
  if (comparison.k_mismatch) {
    log(`(Re-run at k=${comparison.baseline_k} for a valid comparison — see skills_comparison_*.json's k_mismatch field.)`);
  }
}

// ── Main ──────────────────────────────────────────────────────────────────────

await preflightChecks();

log("═".repeat(60));
log(`Skills Progressive-Disclosure Benchmark  run_id=${run_id}`);
log(`k=${CONFIG.k}  condition=simulated_progressive_disclosure_bounded_exploration_forced_synthesis`);
log("═".repeat(60));

log("\nBuilding R function registry for grounding rate...");
const R_FUNCTION_REGISTRY = getRFunctionRegistry();

const recipeDoc = loadRecipeDocument();

let prompts = loadJSONL(CONFIG.prompts_file);
if (CONFIG.prompts !== "all") {
  const wantedIds = CONFIG.prompts.split(",").map(s => s.trim());
  const allIds = new Set(prompts.map(p => p.id));
  const unknownIds = wantedIds.filter(id => !allIds.has(id));
  if (unknownIds.length) {
    log(`❌ Unknown prompt id(s) in --prompts: ${unknownIds.join(", ")}`);
    process.exit(1);
  }
  prompts = prompts.filter(p => wantedIds.includes(p.id));
  log(`  Filtered to ${prompts.length} prompt(s): ${prompts.map(p => p.id).join(", ")}`);
}

// A short TEASER is attached as system_context — not the full recipe text.
// p.prompt is left COMPLETELY unchanged from the vanilla baseline's prompt
// string. Grading uses the same p.id / p.validation / p.expected_workflow
// as the frozen baseline run. The model must discover and read recipe.md
// itself via the view tool (see skills_shim.js) — that's the whole point of
// this script versus recipe_benchmark_runner.js.
const skillTeaser = buildSkillTeaser(recipeDoc.raw);
const skillPrompts = prompts.map(p => ({
  ...p,
  system_context: skillTeaser,
  // p.prompt intentionally untouched
}));

log(`  [skill] teaser size: ${textStats(skillTeaser).chars} chars, ~${textStats(skillTeaser).estimated_tokens} tokens (sent identically to every prompt, every model — recipe.md's ${recipeDoc.stats.chars} chars are NOT pre-sent; the model must view() for them)`);

const modelFilter = CONFIG.model === "all" ? null : CONFIG.model.split(",").map(m => m.trim());
const allResults = {};

for (const [name, fn] of Object.entries(API_MODELS)) {
  if (modelFilter && !modelFilter.includes(name)) continue;
  log(`\n── ${name} (teaser + view tool) ──`);
  const passK   = await runPassK(name, fn, skillPrompts, CONFIG.k);
  const quality = runCodeQuality(name, passK.results, R_FUNCTION_REGISTRY);
  allResults[name] = { passK, quality };
}

// ── Summary ───────────────────────────────────────────────────────────────────

// Aggregate tool_calls_made across every attempt, per model — diagnostic
// signal for which models actually exercise on-demand file access (and for
// what) versus which ones never call the tool even when it's available.
function summarizeToolUsage(passKResults) {
  const pathCounts = {};
  let totalCalls = 0;
  let errorCalls = 0;
  let attemptsWithAnyCall = 0;
  let budgetExhaustedAttempts = 0;
  let totalAttempts = 0;
  for (const promptResult of passKResults) {
    for (const attempt of promptResult.attempts) {
      totalAttempts++;
      const calls = attempt.tool_calls_made || [];
      if (calls.length > 0) attemptsWithAnyCall++;
      if (attempt.tool_budget_exhausted) budgetExhaustedAttempts++;
      totalCalls += calls.length;
      for (const call of calls) {
        // skills_shim.js reports {path, result: "ok"|"error"} objects (view()
        // can fail — nonexistent path, escaped sandbox, binary file — and
        // that's worth tracking, unlike recipe_benchmark_runner.js's narrower
        // tool which could only ever hit a small known allowlist).
        const p = call?.path ?? String(call);
        pathCounts[p] = (pathCounts[p] || 0) + 1;
        if (call?.result === "error") errorCalls++;
      }
    }
  }
  return {
    total_tool_calls: totalCalls,
    error_tool_calls: errorCalls,
    budget_exhausted_attempts: budgetExhaustedAttempts, // needed forced final synthesis — exploration ran up against MAX_TOOL_TURNS
    attempts_using_tool: attemptsWithAnyCall,
    total_attempts: totalAttempts,
    paths_requested: pathCounts, // e.g. { ".": 8, "recipe.md": 8, "code/05_doublets.R": 5, ... }
  };
}

for (const [name, res] of Object.entries(allResults)) {
  if (res.passK) res.toolUsage = summarizeToolUsage(res.passK.results);
}

log("\n" + "─".repeat(68));
log("SKILLS (PROGRESSIVE DISCLOSURE) BENCHMARK SUMMARY");
log("─".repeat(68));

if (Object.keys(allResults).length) {
  log(`${"Model".padEnd(42)} Pass@k     AttemptSR   Workflow    Grounding   ToolUse  Forced`);
  log("─".repeat(104));
  for (const [name, res] of Object.entries(allResults)) {
    if (!res.passK) continue;
    const label = (MODEL_DISPLAY_NAMES[name] || name).padEnd(42);
    const p  = pct(res.passK.avg);
    const asr = pct(res.passK.avgAttemptSuccessRate);
    const w  = pct(res.quality?.avgWorkflow);
    const g  = pct(res.quality?.avgGrounding);
    const tu = res.toolUsage ? `${res.toolUsage.attempts_using_tool}/${res.toolUsage.total_attempts}` : "n/a";
    const fx = res.toolUsage ? String(res.toolUsage.budget_exhausted_attempts) : "n/a";
    log(`${label} ${p.padEnd(10)} ${asr.padEnd(11)} ${w.padEnd(12)} ${g.padEnd(11)} ${tu.padEnd(8)} ${fx}`);
  }
  log("(ToolUse = attempts that made at least 1 view() call / total attempts. Forced = attempts that hit MAX_TOOL_TURNS and got a forced tools-disabled final turn.)");
}

const comparison = buildComparison(allResults);
logComparisonTable(comparison);

// ── Save results ──────────────────────────────────────────────────────────────

const outJson = path.join(CONFIG.results_dir, `skills_results_${run_id}.json`);
const outCsv  = path.join(CONFIG.results_dir, `skills_results_${run_id}_summary.csv`);
const outCompareJson = path.join(CONFIG.results_dir, `skills_comparison_${run_id}.json`);

fs.writeFileSync(outJson, JSON.stringify({
  run_id,
  condition: "simulated_progressive_disclosure_bounded_exploration_forced_synthesis",
  injection_channel: "system_context_teaser_plus_general_view_tool", // metadata teaser only via system context; recipe.md AND supporting files both discovered/read on demand via view() — see skills_shim.js header
  skill_provenance: {
    recipe_file: CONFIG.recipe_file,
    // Hash/stats of the FULL recipe.md — saved for provenance/reproducibility
    // even though it's never pre-sent to the model in this condition. If the
    // model actually called view("recipe.md") during a run, comparing this
    // hash against what the file looked like at run time proves it was
    // reading the same content recipe_benchmark_runner.js's provenance
    // tracks — same file, different access pattern.
    raw_sha256: recipeDoc.sha256,
    stats: recipeDoc.stats,
    // Hash/stats of the SHORT teaser actually sent as systemContext — this
    // is the literal string every model received up front, before any tool
    // calls. Guarantees exact reproducibility of what was visible pre-tool-use.
    teaser_sha256: sha256(skillTeaser),
    teaser_stats: textStats(skillTeaser),
  },
  config: { ...CONFIG },
  results: allResults,
}, null, 2));

function csvCell(value) {
  const s = String(value ?? "");
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

const csvLines = ["model,display_name,pass_k,attempt_success_rate,workflow,grounding"];
for (const [name, res] of Object.entries(allResults)) {
  if (!res.passK) continue;
  const fmt = v => (v === null || v === undefined) ? "" : v.toFixed(4);
  csvLines.push([
    csvCell(name),
    csvCell(MODEL_DISPLAY_NAMES[name] || name),
    fmt(res.passK.avg),
    fmt(res.passK.avgAttemptSuccessRate),
    fmt(res.quality?.avgWorkflow),
    fmt(res.quality?.avgGrounding),
  ].join(","));
}
fs.writeFileSync(outCsv, csvLines.join("\n"));

if (comparison) {
  fs.writeFileSync(outCompareJson, JSON.stringify(comparison, null, 2));
  log(`\nComparison → ${outCompareJson}`);
}

log(`\nResults saved → ${outJson}`);
log(`CSV → ${outCsv}`);
log(`\n✅ Done.`);
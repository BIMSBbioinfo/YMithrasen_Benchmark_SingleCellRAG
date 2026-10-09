// skills_benchmark_runner.js
// "With Skills" condition for the API models. Same prompts, validation and scoring
// as benchmark_runner.js, but the model also gets a short skill teaser in the
// system context and a view(path) tool over recipes/ (tool loop: scripts/skills_shim.js).
//
// The skill is oracle-selected: the teaser says it is relevant, the model still has
// to find and read recipe.md and any supporting files itself. Exploration is capped
// at MAX_TOOL_TURNS; after that the model gets one final turn without tools
// (recorded per attempt as tool_budget_exhausted). This cap is a choice of this
// benchmark, not how Claude Skills work.
//
// Usage (from the repo root):
//   node api_benchmark/skills_benchmark_runner.js --k 3
//   node api_benchmark/skills_benchmark_runner.js --model claude-sonnet-5 --prompts simple_01 --k 1 --skip-comparison
//
// Options:
//   --k <n>                   attempts per prompt (default 3)
//   --model <keys|all>        comma-separated model keys
//   --prompts <ids|all>       comma-separated prompt ids
//   --allow-partial-prompts   run even if prompts lack validation/expected_workflow
//   --baseline <file>         vanilla results to compare against
//   --skip-comparison         do not write the comparison file
//
// Output: skills_results_<run_id>.json, a summary CSV and skills_comparison_<run_id>.json
// in api_benchmark/results/

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

// ---- Config ----

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
    path.resolve(__dirname, "..", "data", "pbmc3k", "filtered_gene_bc_matrices", "hg19")
  ).replace(/\\/g, "/"),

  results_dir:  path.join(__dirname, "results"),
  timeout_ms:   600000, // 10 min

  prompts_file: path.join(__dirname, "prompts/benchmark_prompts.jsonl"),

  // recipes/ sits next to api_benchmark/
  recipe_file:  path.resolve(__dirname, "..", "recipes", "recipe.md"),
  recipe_warning_chars: Number(process.env.BENCH_RECIPE_WARNING_CHARS || 16000),
  recipe_max_chars:     Number(process.env.BENCH_RECIPE_MAX_CHARS     || 20000),

  // vanilla results this run is compared against
  baseline_file: argVal("--baseline",
    path.resolve(__dirname, "..", "results", "api", "vanilla_k3_merged.json")),
};

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 23);
fs.mkdirSync(CONFIG.results_dir, { recursive: true });

// ---- Logging ----

function log(msg) {
  console.log(`[${new Date().toISOString()}] ${msg}`);
}

// ---- Load prompts ----

function loadJSONL(file) {
  return fs.readFileSync(file, "utf8")
    .split("\n").filter(Boolean)
    .map(l => JSON.parse(l));
}

// ---- Call skills shim ----
// The prompt is sent unchanged; the teaser goes separately as systemContext.
// The payload is passed through a temp file because of the Windows command-line length limit.
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
      // up to MAX_TOOL_TURNS API round trips per attempt, so 300 s
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

// ---- Recipe loading & injection ----

function compactMarkdown(text) {
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
    chars: value.length,
    estimated_tokens: estimateTokens(value),
    lines: value.length === 0 ? 0 : value.split(/\r?\n/).length,
  };
}

// Loaded once per run. The raw text is hashed for provenance and used to build the
// teaser; it is never sent to the model directly.
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

// reads a field from recipe.md's YAML frontmatter
function extractFrontmatterField(raw, field) {
  const frontmatterMatch = raw.match(/^---\n([\s\S]*?)\n---/);
  if (!frontmatterMatch) return null;
  const line = frontmatterMatch[1].split("\n").find(l => l.trim().startsWith(`${field}:`));
  return line ? line.split(":").slice(1).join(":").trim() : null;
}

// body of a "## Heading" section, up to the next "##"
function extractSection(raw, heading) {
  const re = new RegExp(`##\\s+${heading.replace(/[.*+?^${}()|[\\]\\\\]/g, "\\$&")}\\s*\\n([\\s\\S]*?)(?=\\n##\\s|$)`, "i");
  const match = raw.match(re);
  return match ? match[1].trim() : null;
}

// Teaser sent as systemContext: skill name and what it does, plus the instruction to
// explore with view(). It does not name recipe.md.
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

// ---- Task-specific validation (identical to benchmark_runner.js) ----

// Turns a prompt's `validation` block into R assertions (same as benchmark_runner.js).
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
    // No object name given: score every Seurat object against the requirements and
    // validate the best match (v1 took the first one by name).
    checks.push(`
.bench_seurat_candidates <- Filter(function(x) {
  v <- tryCatch(get(x, envir = .GlobalEnv), error = function(e) NULL)
  inherits(v, "Seurat")
}, ls(envir = .GlobalEnv))
if (length(.bench_seurat_candidates) == 0) stop("VALIDATION_FAILED: no Seurat object found in environment")
.bench_score_candidate <- function(.obj) {
  .s <- 0
  ${validation.min_cells != null ? `if (ncol(.obj) >= ${validation.min_cells}) .s <- .s + 1` : ""}
  ${validation.min_features != null ? `if (nrow(.obj) >= ${validation.min_features}) .s <- .s + 1` : ""}
  ${(validation.required_reductions || []).map(r => `if ("${r}" %in% names(.obj@reductions)) .s <- .s + 1`).join("\n  ")}
  ${validation.require_clusters ? `if (length(unique(Idents(.obj))) >= ${validation.min_clusters || 2}) .s <- .s + 1` : ""}
  ${(validation.required_columns || []).map(c => `if ("${c}" %in% colnames(.obj@meta.data)) .s <- .s + 1`).join("\n  ")}
  .s
}
.bench_scores <- sapply(.bench_seurat_candidates, function(.nm) {
  .obj <- tryCatch(get(.nm, envir = .GlobalEnv), error = function(e) NULL)
  if (is.null(.obj)) return(-1)
  tryCatch(.bench_score_candidate(.obj), error = function(e) -1)
})
.bench_selected_index <- which.max(.bench_scores)
.bench_obj_name <- .bench_seurat_candidates[[.bench_selected_index]]
.bench_obj <- get(.bench_obj_name, envir = .GlobalEnv)
cat("BENCH_SELECTED_OBJECT:", .bench_obj_name, "\\n")
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

// ---- Execute R code (same as benchmark_runner.js) ----

const PACKAGE_CHECK_PREFIX = `
required_packages <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr","readr","ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(paste0("BENCH_ENV_ERROR: missing required packages: ", paste(missing_packages, collapse = ", "),
              ". Run setup_env.R first."))
}
`.trim();

// Wraps the generated code, runs it with Rscript and reports whether it passed.
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
        // which Seurat object was validated (audit trail)
        const selectedObjectMatch = (stdout || "").match(/BENCH_SELECTED_OBJECT:\s*(\S+)/);
        const selectedObject = selectedObjectMatch ? selectedObjectMatch[1] : null;

        if (!success) {
          log(`  [R exit] error=${error?.code} duration=${duration_ms}ms failureType=${failureType}${selectedObject ? ` selectedObject=${selectedObject}` : ""}`);
          log(`  [R stderr] ${tail(stderr, 600)}`);
          log(`  [R stdout] ${tail(stdout, 400)}`);
        } else {
          if (selectedObject) log(`  [R exit] success selectedObject=${selectedObject}`);
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        resolve({
          success,
          executionSuccess,
          validationSuccess,
          validated: hasValidation,
          failureType,
          selectedObject,
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

// ---- Metrics (same as benchmark_runner.js) ----

const R_REGISTRY_PACKAGES = [
  "base","stats","utils","graphics","grDevices","methods",
  "Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr",
  "readr","ggrepel","DoubletFinder","scDblFinder","SingleCellExperiment"
];

// Returns the set of exported function names, or null if R could not be queried.
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

// Share of function calls that are in the registry or defined in the script.
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

// Workflow score = 0.7 * share of expected steps present + 0.3 * share in the right order.
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

// HumanEval Pass@k estimator; with n = k it is 1 if any attempt passed, else 0.
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

// Fallback extraction from a raw reply: R fence, then any fence, then the whole text.
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

// ---- Pass@k suite ----

// k attempts per prompt for one model. Same grading as benchmark_runner.js; each
// attempt also records the view() calls it made.
async function runPassK(modelName, callerFn, prompts, k) {
  log(`\n[Pass@k] model=${modelName}  k=${k}  prompts=${prompts.length}  condition=skills_progressive_disclosure (teaser + view tool)`);
  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < k; i++) {
      // p.prompt is the unmodified task text; the teaser travels in p.system_context
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
        // separate failure type when the tool budget ran out and the forced final turn
        // still produced no code
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
        selected_object:   exec.selectedObject, // audit trail for object selection
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

// Workflow and grounding for every attempt from runPassK; no extra model calls.
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

// ---- Model registry ----
// No ROLV models here: they do their own retrieval, and skills_shim.js has no ROLV route.

const API_MODELS = {
  "gpt-5.5":           (p, sc) => callBaseline("gpt-5.5", p, sc),
  "gpt-5.4-mini":      (p, sc) => callBaseline("gpt-5.4-mini", p, sc),
  "gemini-openrouter": (p, sc) => callBaseline("gemini-openrouter", p, sc),
  "gemini-flash-lite": (p, sc) => callBaseline("gemini-flash-lite", p, sc),
  "claude-opus-4-8":   (p, sc) => callBaseline("claude-opus-4-8", p, sc),
  "claude-sonnet-5":   (p, sc) => callBaseline("claude-sonnet-5", p, sc),
  "claude-haiku-4-5":  (p, sc) => callBaseline("claude-haiku-4-5", p, sc),
};

const MODEL_DISPLAY_NAMES = {
  "gpt-5.5":           "GPT-5.5 (OpenAI flagship) + skill (view tool)",
  "gpt-5.4-mini":      "GPT-5.4 Mini (OpenAI) + skill (view tool)",
  "gemini-openrouter": "Gemini 3.5 Flash (via OpenRouter) + skill (view tool)",
  "gemini-flash-lite": "Gemini 3.1 Flash-Lite (via OpenRouter) + skill (view tool)",
  "claude-opus-4-8":   "Claude Opus 4.8 (Anthropic direct) + skill (view tool)",
  "claude-sonnet-5":   "Claude Sonnet 5 (Anthropic direct) + skill (view tool)",
  "claude-haiku-4-5":  "Claude Haiku 4.5 (Anthropic direct) + skill (view tool)",
};

const MODEL_REQUIRED_KEYS = {
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
  ANTHROPIC_API_KEY:  "Anthropic key",
};

// ---- Preflight checks ----

// Checks keys, recipe, R, packages, dataset and prompts, then tests the teaser and
// the view tool with one live call per model. Exits if something is wrong.
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

  // baseline file for the comparison (warning only)
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
        // loading Seurat in a fresh Rscript is slow, hence 60 s
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
    // The smoke test checks three things per model: the teaser reaches the model, view()
    // works for both listing and reading, and the tool result shows up in the answer.
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

// ---- Comparison against frozen baseline ----

function pct(v) { return v === null || v === undefined ? "N/A" : (v * 100).toFixed(1) + "%"; }

// Per-model and per-prompt difference to the vanilla baseline file.
function buildComparison(skillResults) {
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

  // Pass@k is not comparable across different k, so a mismatch is flagged, not hidden.
  const baselineK = baseline.config?.k;
  const kMismatch = baselineK != null && baselineK !== CONFIG.k;
  if (kMismatch) {
    log(`\n⚠️  K MISMATCH: this run used k=${CONFIG.k}, baseline used k=${baselineK}. Pass@k deltas below are NOT a like-for-like comparison — Pass@${baselineK} is structurally >= Pass@${CONFIG.k} for the same success rate, independent of any skill effect. Treat deltas as directional/sanity-check only until re-run at matching k.`);
  }

  const rows = [];
  for (const [modelName, skillRes] of Object.entries(skillResults)) {
    const baseRes = baseline.results?.[modelName];
    if (!baseRes?.passK) {
      log(`  [compare] no baseline entry for model="${modelName}" — skipping`);
      continue;
    }
    const basePassK   = baseRes.passK.avg;
    const baseASR      = baseRes.passK.avgAttemptSuccessRate;
    const baseWorkflow = baseRes.quality?.avgWorkflow ?? null;
    const baseGrounding= baseRes.quality?.avgGrounding ?? null;

    const skillPassK   = skillRes.passK.avg;
    const skillASR      = skillRes.passK.avgAttemptSuccessRate;
    const skillWorkflow = skillRes.quality?.avgWorkflow ?? null;
    const skillGrounding= skillRes.quality?.avgGrounding ?? null;

    rows.push({
      model: modelName,
      baseline_pass_k: basePassK,
      skill_pass_k: skillPassK,
      delta_pass_k: skillPassK - basePassK,
      baseline_attempt_success_rate: baseASR,
      skill_attempt_success_rate: skillASR,
      delta_attempt_success_rate: skillASR - baseASR,
      baseline_workflow: baseWorkflow,
      skill_workflow: skillWorkflow,
      delta_workflow: (baseWorkflow != null && skillWorkflow != null) ? skillWorkflow - baseWorkflow : null,
      baseline_grounding: baseGrounding,
      skill_grounding: skillGrounding,
      delta_grounding: (baseGrounding != null && skillGrounding != null) ? skillGrounding - baseGrounding : null,
      // per-prompt breakdown
      per_prompt: (skillRes.passK.results || []).map(sp => {
        const bp = (baseRes.passK.results || []).find(x => x.id === sp.id);
        return {
          id: sp.id,
          baseline_score: bp ? bp.score : null,
          skill_score: sp.score,
          delta: bp ? sp.score - bp.score : null,
        };
      }),
    });
  }

  return {
    baseline_run_id: baseline.run_id,
    baseline_file: CONFIG.baseline_file,
    skill_k: CONFIG.k,
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
    log(`⚠️  K MISMATCH: this run k=${comparison.skill_k} vs baseline k=${comparison.baseline_k} — deltas below are NOT a valid Pass@k comparison, directional only.`);
  }
  log("─".repeat(68));
  log(`${"Model".padEnd(30)} Baseline    +Skill      Δ Pass@k`);
  log("─".repeat(68));
  for (const r of comparison.rows) {
    const label = (MODEL_DISPLAY_NAMES[r.model] || r.model).padEnd(30);
    const deltaStr = (r.delta_pass_k >= 0 ? "+" : "") + (r.delta_pass_k * 100).toFixed(1) + "pp";
    log(`${label} ${pct(r.baseline_pass_k).padEnd(11)} ${pct(r.skill_pass_k).padEnd(11)} ${deltaStr}`);
  }
  log("─".repeat(68));
  log("(pp = percentage points. Positive Δ = the skill-enabled condition helped this model on this metric.)");
  if (comparison.k_mismatch) {
    log(`(Re-run at k=${comparison.baseline_k} for a valid comparison — see skills_comparison_*.json's k_mismatch field.)`);
  }
}

// ---- Main ----

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

// attach the teaser as system_context; p.prompt stays as in the vanilla run
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

// ---- Summary ----

// view() calls per model: how many, which paths, how many attempts hit the turn limit
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
        // the shim reports {path, result: "ok"|"error"} per call
        const p = call?.path ?? String(call);
        pathCounts[p] = (pathCounts[p] || 0) + 1;
        if (call?.result === "error") errorCalls++;
      }
    }
  }
  return {
    total_tool_calls: totalCalls,
    error_tool_calls: errorCalls,
    budget_exhausted_attempts: budgetExhaustedAttempts, // hit MAX_TOOL_TURNS
    attempts_using_tool: attemptsWithAnyCall,
    total_attempts: totalAttempts,
    paths_requested: pathCounts, // e.g. { ".": 8, "recipe.md": 8, "code/05_doublets.R": 5 }
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

// ---- Save results ----

const outJson = path.join(CONFIG.results_dir, `skills_results_${run_id}.json`);
const outCsv  = path.join(CONFIG.results_dir, `skills_results_${run_id}_summary.csv`);
const outCompareJson = path.join(CONFIG.results_dir, `skills_comparison_${run_id}.json`);

fs.writeFileSync(outJson, JSON.stringify({
  run_id,
  // Bump when buildValidationR changes pass/fail outcomes; results are not comparable
  // across versions.
  grader_version: "validation_object_selection_v2",
  condition: "simulated_progressive_disclosure_bounded_exploration_forced_synthesis",
  injection_channel: "system_context_teaser_plus_general_view_tool", // teaser in the system context, files read on demand via view()
  skill_provenance: {
    recipe_file: CONFIG.recipe_file,
    // hash of the full recipe.md, for provenance (it is never pre-sent to the model)
    raw_sha256: recipeDoc.sha256,
    stats: recipeDoc.stats,
    // hash of the teaser every model received
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
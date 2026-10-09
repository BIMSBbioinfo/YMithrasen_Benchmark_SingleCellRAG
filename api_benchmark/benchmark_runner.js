// benchmark_runner.js
// Runs the API models without Skills and the Agentic RAG (ROLV) systems on the
// benchmark prompts, executes the generated R code and scores it.
//
// Usage (from the repo root):
//   node api_benchmark/benchmark_runner.js --k 3 --model gpt-5.5,claude-sonnet-5
//
// Options:
//   --k <n>                   attempts per prompt (default 3)
//   --model <keys|all>        comma-separated model keys
//   --prompts <ids|all>       comma-separated prompt ids
//   --allow-partial-prompts   run even if prompts lack validation/expected_workflow
//
// Output: results_<run_id>.json and results_<run_id>_summary.csv in api_benchmark/results/

import "dotenv/config";
import fs   from "fs";
import path from "path";
import { execFile, exec, execFileSync } from "child_process";
import { fileURLToPath }  from "url";
import { createReadStream } from "fs";
import readline from "readline";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ---- Config ----

const CONFIG = {
  k:            parseInt(process.argv.find((a, i) => process.argv[i-1] === "--k") || "3"),
  suite:        process.argv.find((a, i) => process.argv[i-1] === "--suite") || "internal",
  model:        process.argv.find((a, i) => process.argv[i-1] === "--model") || "all",
  prompts:      process.argv.find((a, i) => process.argv[i-1] === "--prompts") || "all",
  allowPartialPrompts: process.argv.includes("--allow-partial-prompts"),
  r_exec:       process.env.R_EXEC || "Rscript", // Rscript from PATH unless R_EXEC is set
  dataset_dir:  (process.env.BENCH_DATASET_DIR ||
    path.resolve(__dirname, "..", "data", "pbmc3k", "filtered_gene_bc_matrices", "hg19")
  ).replace(/\\/g, "/"),
  results_dir:  path.join(__dirname, "results"),
  timeout_ms:   600000, // 10 min
  prompts_file: path.join(__dirname, "prompts/benchmark_prompts.jsonl"),
};

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 23);
fs.mkdirSync(CONFIG.results_dir, { recursive: true });

// ---- Logging ----

function log(msg) {
  const line = `[${new Date().toISOString()}] ${msg}`;
  console.log(line);
}

// ---- Load prompts ----

function loadJSONL(file) {
  return fs.readFileSync(file, "utf8")
    .split("\n").filter(Boolean)
    .map(l => JSON.parse(l));
}

// ---- Call baseline shim ----

// Calls one model through baseline_shim.js in a child process.
// Never rejects; failures come back in `error`.
function callBaseline(model, prompt) {
  return new Promise((resolve) => {
    const payload = Buffer.from(JSON.stringify({ model, prompt })).toString("base64");
    const shimPath = path.join(__dirname, "scripts/baseline_shim.js");
    exec(
      `node "${shimPath}" "${payload}"`,
      { timeout: 180000 }, // reasoning models need more than 60 s
      (error, stdout) => {
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

// ---- Call ROLV orchestrator shim ----

// Calls the ROLV orchestrator. shimFile selects the variant with or without the
// pre-execution planner.
function callRolv(prompt, shimFile = "orchestrator_shim.js") {
  return new Promise((resolve) => {
    const shimPath = path.join(__dirname, "scripts", shimFile);
    const payload  = Buffer.from(JSON.stringify({ prompt })).toString("base64");
    exec(
      `node "${shimPath}" "${payload}"`,
      { timeout: 90000 },
      (error, stdout, stderr) => {
        try {
          const lines  = (stdout || "").trim().split("\n").filter(l => l.startsWith("{"));
          const parsed = JSON.parse(lines.pop() || "{}");
          resolve({ result: parsed.result || parsed.response || null, error: error?.message || null });
        } catch {
          resolve({ result: null, error: error?.message || "parse error" });
        }
      }
    );
  });
}

// ---- Task-specific validation ----
// Turns a prompt's `validation` block into R assertions that run after the
// generated code. Prompts without one are scored on execution only.
function buildValidationR(validation = {}) {
  const checks = [];
  const explicitObj = validation.primary_object;
  // only these fields need the primary object
  const needsPrimaryObject =
    validation.min_cells != null ||
    validation.min_features != null ||
    validation.required_reductions?.length ||
    validation.require_clusters ||
    validation.required_columns?.length;

  // No object name given: score every Seurat object in the environment against the
  // requirements and validate the best match. (v1 took the first one by name, which
  // picked a raw `seu` over a processed `seu_filtered`.)
  if (needsPrimaryObject && !explicitObj) {
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

// ---- Execute R code ----

// Packages are installed once with the environment (environment.yml, setup_env.R),
// never during a run. Missing packages fail fast.
const PACKAGE_CHECK_PREFIX = `
required_packages <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr","readr","ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(paste0("BENCH_ENV_ERROR: missing required packages: ", paste(missing_packages, collapse = ", "),
              ". Run setup_env.R first."))
}
`.trim();

// Wraps the generated code (setwd, package check, tryCatch, validation), runs it
// with Rscript and reports whether it passed. The temp dir is kept on failure.
function executeRCode(rCode, datasetDir, validation = {}) {
  const tmpId   = Math.random().toString(36).slice(2, 8);
  const tmpDir  = path.join(CONFIG.results_dir, `rexec_${tmpId}`);
  fs.mkdirSync(tmpDir, { recursive: true });
  const scriptPath = path.join(tmpDir, "bench.R");

  // Normalise all placeholder and tutorial paths to "."
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
        // no validation defined: success means execution only
        const validationSuccess = hasValidation ? (stdout || "").includes("BENCH_VALIDATION_OK") : true;
        const success = executionSuccess && validationSuccess;
        const failureType = success ? null : classifyRFailure(stderr, stdout, error);
        // which Seurat object was validated (audit trail)
        const selectedObjectMatch = (stdout || "").match(/BENCH_SELECTED_OBJECT:\s*(\S+)/);
        const selectedObject = selectedObjectMatch ? selectedObjectMatch[1] : null;

        if (!success) {
          // log the tail: R errors come after Seurat's progress output
          log(`  [R exit] error=${error?.code} duration=${duration_ms}ms failureType=${failureType}${selectedObject ? ` selectedObject=${selectedObject}` : ""}`);
          log(`  [R stderr] ${tail(stderr, 600)}`);
          log(`  [R stdout] ${tail(stdout, 400)}`);
        } else {
          if (selectedObject) log(`  [R exit] success selectedObject=${selectedObject}`);
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        // the full stderr goes into the results JSON
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

// Separates environment problems from errors in the generated code.
// Order matters: specific patterns first.
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

// ---- Metrics ----
// Grounding is checked against the functions exported by the installed R packages,
// queried from R itself.

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

  // run from a file: a long script passed via -e gets mangled on Windows
  const tmpPath = path.join(CONFIG.results_dir, `_registry_query_${Date.now()}.R`);
  fs.mkdirSync(CONFIG.results_dir, { recursive: true });
  fs.writeFileSync(tmpPath, registryScript);

  try {
    const output = execFileSync(CONFIG.r_exec, ["--vanilla", tmpPath],
      { encoding: "utf8", timeout: 120000 });
    const entries = output.split("\n").map(l => l.trim()).filter(Boolean);
    const set = new Set(entries);
    log(`  R function registry built: ${set.size} exported functions from ${R_REGISTRY_PACKAGES.length} packages`);
    if (set.size < 100) {
      log(`  ⚠️  Registry looks suspiciously small (${set.size}) — grounding rate for this run may be unreliable. Sample: ${entries.slice(0, 10).join(", ") || "(empty)"}`);
    } else {
      log(`  Sample entries: ${entries.slice(0, 8).join(", ")}...`);
    }
    return set;
  } catch (e) {
    log(`  ⚠️  Could not build R function registry: ${e.message} — grounding rate will be unavailable this run`);
    return null; // null = registry unavailable, distinct from empty set
  } finally {
    try { fs.unlinkSync(tmpPath); } catch {}
  }
}

const NON_FUNCTION_KEYWORDS = new Set([
  "if","else","for","while","repeat","function","return","switch","next","break"
]);

// functions the script defines itself, so local helpers do not count as ungrounded
function extractDefinedFunctions(code) {
  const defined = new Set();
  const pattern = /\b([A-Za-z][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(/g;
  for (const m of code.matchAll(pattern)) defined.add(m[1]);
  return defined;
}

// Share of function calls that are in the registry or defined in the script.
// null if the registry is unavailable.
function calcGroundingRate(code, registry) {
  if (!code) return { rate: 0, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: !!registry };
  if (!registry) return { rate: null, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: false };

  code = stripRComments(code);
  const definedFunctions = extractDefinedFunctions(code);
  const calls = [...code.matchAll(/\b([A-Za-z][A-Za-z0-9_.]*)\s*\(/g)]
    .map(m => m[1])
    .filter(name => !NON_FUNCTION_KEYWORDS.has(name));

  if (!calls.length) return { rate: 1, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: true };

  const grounded   = [];
  const ungrounded = [];
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

function stripRComments(code) {
  // does not handle # inside strings (rare in generated Seurat code)
  return code.split("\n").map(line => line.replace(/#.*$/, "")).join("\n");
}

// Workflow score = 0.7 * share of expected steps present + 0.3 * share in the right order.
// null if the prompt has no expected_workflow.
function calcWorkflowScore(code, expectedWorkflow = []) {
  if (!expectedWorkflow.length) return { score: null, presenceScore: null, orderScore: null, steps: [] };
  if (!code) return { score: 0, presenceScore: 0, orderScore: 0, steps: [] };

  code = stripRComments(code);

  const evaluated = expectedWorkflow.map(stepName => {
    const step = STEP_MAP.get(stepName);
    if (!step) {
      // unknown step name (typo): flagged so it is not read as a model failure
      return { name: stepName, present: false, position: -1, unknown: true };
    }
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

// HumanEval Pass@k estimator, 1 - C(n-c, k) / C(n, k).
// With n = k this is 1 if any attempt passed, else 0.
function passAtK(n, c, k) {
  if (n - c < k) return 1.0;
  let num = 1, den = 1;
  for (let i = 0; i < k; i++) {
    num *= (n - c - i);
    den *= (n - i);
  }
  return 1 - num / den;
}

// used to check that the k attempts actually differ
function normalizeCodeForComparison(code) {
  return (code || "").replace(/#.*$/gm, "").replace(/\s+/g, " ").trim();
}

function countUniqueGenerations(attempts) {
  const unique = new Set(attempts.map(a => normalizeCodeForComparison(a.code)).filter(Boolean));
  return unique.size;
}

// ---- Pass@k suite ----

// k attempts per prompt for one model: generate, execute, validate.
async function runPassK(modelName, callerFn, prompts, k) {
  log(`\n[Pass@k] model=${modelName}  k=${k}  prompts=${prompts.length}`);
  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < k; i++) {
      const res  = await callerFn(p.prompt);
      // the shim already extracted the code; extracting again rejected valid short snippets
      const code = res.alreadyExtracted
        ? (res.result || null)
        : extractRFromText(res.result || "");
      if (!code) {
        log(`  ⚠️  [${p.id}] no_r_code_extracted — model_error=${res.error || "(none)"} raw_preview=${JSON.stringify((res.raw_preview || "(not captured)").slice(0, 300))}`);
        attempts.push({ passed: false, error: "no_r_code_extracted", failure_type: "no_code_extracted", model_error: res.error || null, raw_preview: res.raw_preview || null });
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

// Fallback extraction from a raw reply: R fence, then any fence, then the whole
// text if it looks like R.
function extractRFromText(text) {
  if (!text) return null;
  // Strip ROLV scratchpad blocks
  text = text.replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "").trim();
  // Try R-specific fence first
  let match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  // Fall back to any code fence
  match = text.match(/```[a-z]*\n?([\s\S]*?)```/i);
  if (match) return match[1].trim();
  // no fence: take the whole text if it looks like R (rolv-raw returns plain code)
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) {
    return text;
  }
  return null;
}

// ---- Code quality ----
// Workflow and grounding for every attempt from runPassK; no extra model calls.
function runCodeQuality(modelName, passKResults, registry) {
  log(`\n[CodeQuality] model=${modelName} (scored from all Pass@k attempts, registry ${registry ? "available" : "UNAVAILABLE"})`);
  const scores = [];

  for (const r of passKResults) {
    const attemptScores = r.attempts.map((attempt, i) => {
      const code = attempt.code || null;
      const workflowResult  = code ? calcWorkflowScore(code, r.expected_workflow) : { score: 0, presenceScore: 0, orderScore: 0, steps: [] };
      const groundingResult = code ? calcGroundingRate(code, registry)            : { rate: 0, groundedCalls: [], ungroundedCalls: [], totalCalls: 0, registryAvailable: !!registry };
      return { attempt: i + 1, passed: attempt.passed, workflow: workflowResult, grounding: groundingResult };
    });

    // scores can be null (no expected_workflow, or no registry)
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

  log(`  avg workflow success rate: ${avgWorkflow !== null ? avgWorkflow.toFixed(3) : "N/A (no prompts declared expected_workflow)"}`);
  log(`  avg grounding rate:        ${avgGrounding !== null ? avgGrounding.toFixed(3) : "N/A (registry unavailable)"}`);
  return { avgWorkflow, avgGrounding, scores };
}

// ---- Preflight checks ----

// API keys each model needs. Preflight only checks the selected models.
const MODEL_REQUIRED_KEYS = {
  "rolv":                ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "rolv-nopreexec":      ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "rolv-raw":            ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "gpt-5.5":             ["OPENAI_API_KEY"],
  "gpt-5.4-mini":        ["OPENAI_API_KEY"],
  "gemini-openrouter":   ["OPENROUTER_API_KEY"],
  "gemini-flash-lite":   ["OPENROUTER_API_KEY"],
  "claude-opus-4-8":     ["ANTHROPIC_API_KEY"],
  "claude-sonnet-5":     ["ANTHROPIC_API_KEY"],
  "claude-haiku-4-5":    ["ANTHROPIC_API_KEY"],
};

const KEY_LABELS = {
  ROLV_API_URL:       "ROLV API URL",
  ROLV_JWT_TOKEN:     "ROLV JWT token",
  OPENAI_API_KEY:     "OpenAI key",
  OPENROUTER_API_KEY: "OpenRouter key",
  ANTHROPIC_API_KEY:  "Anthropic key",
};

// Checks keys, R, packages, dataset and prompts, then makes one live call per model.
// Exits before any prompt is run if something is wrong.
async function preflightChecks() {
  log("═".repeat(60));
  log("PREFLIGHT CHECKS");
  log("═".repeat(60));

  const errors   = [];
  const warnings = [];

  const selectedModels = CONFIG.model === "all"
    ? Object.keys(MODEL_REQUIRED_KEYS)
    : CONFIG.model.split(",").map(m => m.trim());
  const needsRolv = selectedModels.some(m => ["rolv", "rolv-nopreexec", "rolv-raw"].includes(m));

  log(`  Models in scope: ${selectedModels.join(", ")}`);

  // catch typos in --model before they silently drop a model
  const unknownModels = selectedModels.filter(
    m => !Object.prototype.hasOwnProperty.call(MODEL_REQUIRED_KEYS, m)
  );
  if (unknownModels.length) {
    errors.push(
      `Unknown model key(s): ${unknownModels.join(", ")} — check spelling against: ${Object.keys(MODEL_REQUIRED_KEYS).join(", ")}`
    );
  }

  // 1. API keys for the selected models
  const neededKeys = new Set();
  for (const m of selectedModels) {
    for (const key of (MODEL_REQUIRED_KEYS[m] || [])) neededKeys.add(key);
  }

  for (const key of neededKeys) {
    const val = process.env[key];
    if (!val || val.startsWith("sk-replace") || val.startsWith("sk-or-replace")) {
      errors.push(`Missing or placeholder: ${key} (${KEY_LABELS[key] || key})`);
    } else {
      log(`  ✅ ${key} present`);
    }
  }
  const skippedKeys = Object.values(KEY_LABELS).length
    ? Object.keys(KEY_LABELS).filter(k => !neededKeys.has(k))
    : [];
  if (skippedKeys.length) {
    log(`  ⏭️  Skipping unused keys for this run: ${skippedKeys.join(", ")}`);
  }

  // 2. JWT expiry, only if a ROLV model is selected
  const jwt = process.env.ROLV_JWT_TOKEN;
  if (needsRolv && jwt && !jwt.startsWith("sk-replace")) {
    try {
      const payload = JSON.parse(Buffer.from(jwt.split(".")[1], "base64").toString("utf8"));
      const expiresAt = new Date(payload.exp * 1000);
      if (Date.now() > payload.exp * 1000) {
        errors.push(`ROLV_JWT_TOKEN is EXPIRED (expired ${expiresAt.toISOString()}) — refresh it before running`);
      } else {
        log(`  ✅ ROLV_JWT_TOKEN valid until ${expiresAt.toISOString()}`);
      }
    } catch {
      warnings.push("Could not decode ROLV_JWT_TOKEN — check it is a valid JWT");
    }
  }

  // 3. R and packages. Runs Rscript instead of fs.existsSync because r_exec may be a PATH command.
  try {
    const versionOut = execFileSync(CONFIG.r_exec, ["--version"], { encoding: "utf8", timeout: 15000 });
    log(`  ✅ R reachable: ${versionOut.split("\n")[0].trim()}`);

    const pkgCheckScript = `
required_packages <- c("Seurat","SeuratObject","dplyr","ggplot2","patchwork","Matrix","tidyr","stringr","readr","ggrepel","DoubletFinder","scDblFinder","SingleCellExperiment")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ":", as.character(packageVersion(pkg)), "\\n")
  } else {
    cat(pkg, ": MISSING\\n")
  }
}
`.trim();
    const pkgTmpPath = path.join(CONFIG.results_dir, `_pkgcheck_${Date.now()}.R`);
    fs.mkdirSync(CONFIG.results_dir, { recursive: true });
    fs.writeFileSync(pkgTmpPath, pkgCheckScript);
    let pkgOut;
    try {
      // loading Seurat in a fresh Rscript takes about 35 s, hence 60 s
      pkgOut = execFileSync(CONFIG.r_exec, ["--vanilla", pkgTmpPath], { encoding: "utf8", timeout: 60000 });
    } finally {
      try { fs.unlinkSync(pkgTmpPath); } catch {}
    }
    const missing = [];
    for (const line of pkgOut.split("\n").map(l => l.trim()).filter(Boolean)) {
      log(`    ${line}`);
      if (line.endsWith("MISSING")) missing.push(line.split(":")[0].trim());
    }
    if (missing.length) {
      errors.push(`Missing R packages: ${missing.join(", ")} — run setup_env.R first`);
    }
  } catch (e) {
    errors.push(`R not reachable via "${CONFIG.r_exec}": ${e.message} — install R or set R_EXEC`);
  }

  // 4. Dataset directory exists
  if (!fs.existsSync(CONFIG.dataset_dir)) {
    errors.push(`Dataset directory not found: ${CONFIG.dataset_dir}`);
  } else {
    log(`  ✅ Dataset directory found: ${CONFIG.dataset_dir}`);
  }

  // 5. Prompts file and annotations. Missing validation/expected_workflow is an error
  // unless --allow-partial-prompts is given.
  if (!fs.existsSync(CONFIG.prompts_file)) {
    errors.push(`Prompts file not found: ${CONFIG.prompts_file}`);
  } else {
    log(`  ✅ Prompts file found`);
    if (CONFIG.suite === "internal" || CONFIG.suite === "all") {
      try {
        const promptsPreview = loadJSONL(CONFIG.prompts_file);
        const missingValidation = promptsPreview.filter(p => !p.validation || Object.keys(p.validation).length === 0);
        const missingWorkflow   = promptsPreview.filter(p => !p.expected_workflow || p.expected_workflow.length === 0);
        const bucket = CONFIG.allowPartialPrompts ? warnings : errors;

        if (missingValidation.length) {
          bucket.push(
            `${missingValidation.length}/${promptsPreview.length} prompts have no 'validation' field — ` +
            `Pass@k would be execution-only, not task-validated, for: ${missingValidation.map(p => p.id).join(", ")}` +
            (CONFIG.allowPartialPrompts ? "" : " (pass --allow-partial-prompts to run anyway)")
          );
        }
        if (missingWorkflow.length) {
          bucket.push(
            `${missingWorkflow.length}/${promptsPreview.length} prompts have no 'expected_workflow' — workflow score would be N/A for: ${missingWorkflow.map(p => p.id).join(", ")}` +
            (CONFIG.allowPartialPrompts ? "" : " (pass --allow-partial-prompts to run anyway)")
          );
        }

        // an unknown step name is always an error
        const knownSteps = new Set(PIPELINE_STEPS.map(s => s.name));
        const badSteps = [];
        for (const p of promptsPreview) {
          for (const step of (p.expected_workflow || [])) {
            if (!knownSteps.has(step)) badSteps.push(`${p.id}:"${step}"`);
          }
        }
        if (badSteps.length) {
          errors.push(
            `Unknown expected_workflow step name(s) — check for typos against PIPELINE_STEPS (${[...knownSteps].join(", ")}): ${badSteps.join(", ")}`
          );
        }
      } catch (e) {
        warnings.push(`Could not pre-check prompt validation coverage: ${e.message}`);
      }
    }
  }

  // 6. One live call per selected model
  if (!errors.length) {
    log("\n  Running live API smoke tests (1 call each, selected models only)...");
    const SMOKE_PROMPT = "Write a single line of R code that assigns the number 1 to a variable called x.";

    const ALL_SMOKE_TESTS = {
      "gpt-5.5":            () => callBaseline("gpt-5.5", SMOKE_PROMPT),
      "gpt-5.4-mini":       () => callBaseline("gpt-5.4-mini", SMOKE_PROMPT),
      "gemini-openrouter":  () => callBaseline("gemini-openrouter", SMOKE_PROMPT),
      "gemini-flash-lite":  () => callBaseline("gemini-flash-lite", SMOKE_PROMPT),
      "claude-opus-4-8":    () => callBaseline("claude-opus-4-8", SMOKE_PROMPT),
      "claude-sonnet-5":    () => callBaseline("claude-sonnet-5", SMOKE_PROMPT),
      "claude-haiku-4-5":   () => callBaseline("claude-haiku-4-5", SMOKE_PROMPT),
      "rolv": () => callRolv(SMOKE_PROMPT), // covers rolv, rolv-nopreexec, rolv-raw — same endpoint
    };

    const smokeLabels = selectedModels.map(m => ["rolv", "rolv-nopreexec", "rolv-raw"].includes(m) ? "rolv" : m);
    const uniqueSmokeLabels = [...new Set(smokeLabels)].filter(l => ALL_SMOKE_TESTS[l]);

    for (const label of uniqueSmokeLabels) {
      try {
        const res = await ALL_SMOKE_TESTS[label]();
        const hasResult = res.result != null && res.result !== "";
        if (res.error && !hasResult) throw new Error(res.error);
        log(`  ✅ ${label} reachable${hasResult ? "" : " (no content returned — check API key)"}`);
      } catch (e) {
        errors.push(`${label} smoke test FAILED: ${e.message}`);
        log(`  ❌ ${label}: ${e.message}`);
      }
    }
  }

  // 7. Report
  log("\n" + "─".repeat(60));
  if (errors.length) {
    log("❌ PREFLIGHT FAILED — fix these before running:");
    errors.forEach(e => log(`   • ${e}`));
    log("\nAborting. No prompts were run.");
    process.exit(1);
  }
  if (warnings.length) {
    warnings.forEach(w => log(`  ⚠️  ${w}`));
  }
  log("✅ All preflight checks passed — starting benchmark run.");
  log("─".repeat(60));
}

await preflightChecks();

log("═".repeat(60));
log(`ROLV Benchmark  run_id=${run_id}`);
log(`k=${CONFIG.k}  suite=${CONFIG.suite}`);
log("═".repeat(60));

// built once, used for every model's grounding score
log("\nBuilding R function registry for grounding rate...");
const R_FUNCTION_REGISTRY = getRFunctionRegistry();

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

const API_MODELS = {
  rolv:                (p) => callRolv(p),
  "rolv-nopreexec":    (p) => callRolv(p, "orchestrator_shim_nopreexec.js"),
  "rolv-raw":          (p) => callBaseline("rolv-raw", p),
  "gpt-5.5":           (p) => callBaseline("gpt-5.5", p),
  "gpt-5.4-mini":      (p) => callBaseline("gpt-5.4-mini", p),
  "gemini-openrouter": (p) => callBaseline("gemini-openrouter", p),
  "gemini-flash-lite": (p) => callBaseline("gemini-flash-lite", p),
  "claude-opus-4-8":   (p) => callBaseline("claude-opus-4-8", p),
  "claude-sonnet-5":   (p) => callBaseline("claude-sonnet-5", p),
  "claude-haiku-4-5":  (p) => callBaseline("claude-haiku-4-5", p),
};

const MODEL_DISPLAY_NAMES = {
  "rolv":              "ROLV (Claude Sonnet + RAG + PreExec)",
  "rolv-nopreexec":    "ROLV (Claude Sonnet + RAG only)",
  "rolv-raw":          "Claude Sonnet 4.5 (via ROLV raw, no RAG)",
  "gpt-5.5":           "GPT-5.5 (OpenAI flagship)",
  "gpt-5.4-mini":      "GPT-5.4 Mini (OpenAI)",
  "gemini-openrouter": "Gemini 3.5 Flash (via OpenRouter)",
  "gemini-flash-lite": "Gemini 3.1 Flash-Lite (via OpenRouter)",
  "claude-opus-4-8":   "Claude Opus 4.8 (Anthropic direct)",
  "claude-sonnet-5":   "Claude Sonnet 5 (Anthropic direct)",
  "claude-haiku-4-5":  "Claude Haiku 4.5 (Anthropic direct)",
};

const allResults = {};

if (CONFIG.suite === "internal" || CONFIG.suite === "all") {
  log(`\n══ Internal Suite ══`);
  const modelFilter = CONFIG.model === "all"
    ? null
    : CONFIG.model.split(",").map(m => m.trim());
  for (const [name, fn] of Object.entries(API_MODELS)) {
    if (modelFilter && !modelFilter.includes(name)) continue;
    log(`\n── ${name} ──`);
    const passK   = await runPassK(name, fn, prompts, CONFIG.k);
    const quality = runCodeQuality(name, passK.results, R_FUNCTION_REGISTRY);
    allResults[name] = { passK, quality };
  }
}

// ---- Summary ----

log("\n" + "─".repeat(68));
log("BENCHMARK SUMMARY");
log("─".repeat(68));

function pct(v) {
  return v === null || v === undefined ? "N/A" : (v * 100).toFixed(1) + "%";
}

if (allResults && Object.keys(allResults).length) {
  log("\nINTERNAL METRICS");
  log(`${"Model".padEnd(42)} Pass@k     AttemptSR   Workflow    Grounding`);
  log("─".repeat(84));
  for (const [name, res] of Object.entries(allResults)) {
    if (!res.passK) continue;
    const label = (MODEL_DISPLAY_NAMES[name] || name).padEnd(42);
    const p  = pct(res.passK.avg);
    const asr = pct(res.passK.avgAttemptSuccessRate);
    const w  = pct(res.quality?.avgWorkflow);
    const g  = pct(res.quality?.avgGrounding);
    log(`${label} ${p.padEnd(10)} ${asr.padEnd(11)} ${w.padEnd(12)} ${g}`);
  }
  log("\n(AttemptSR = successful attempts / total attempts — distinguishes 1/3-pass from 3/3-pass models with identical Pass@k)");
}

// Save results
const outJson = path.join(CONFIG.results_dir, `results_${run_id}.json`);
const outCsv  = path.join(CONFIG.results_dir, `results_${run_id}_summary.csv`);
const MODEL_SAMPLING = {
  "rolv":                { temperature: 0.6 },
  "rolv-nopreexec":      { temperature: 0.6 },
  "rolv-raw":            { temperature: 0.6 },
  "gpt-5.5":              { temperature: null, reason: "model does not accept temperature parameter" },
  "gpt-5.4-mini":        { temperature: 0.6 },
  "gemini-openrouter":   { temperature: 0.6 },
  "gemini-flash-lite":   { temperature: 0.6 },
  "claude-opus-4-8":     { temperature: null, reason: "model does not accept temperature parameter" },
  "claude-sonnet-5":     { temperature: null, reason: "model does not accept temperature parameter" },
  "claude-haiku-4-5":    { temperature: null, reason: "model does not accept temperature parameter" },
};

fs.writeFileSync(outJson, JSON.stringify({
  run_id,
  // Bump when buildValidationR changes pass/fail outcomes; results are not comparable
  // across versions. v2 validates the best-matching Seurat object instead of the first by name.
  grader_version: "validation_object_selection_v2",
  config: { ...CONFIG, grounding_registry_packages: R_REGISTRY_PACKAGES, model_sampling: MODEL_SAMPLING },
  results: allResults,
}, null, 2));

// CSV
const csvLines = ["model,display_name,pass_k,attempt_success_rate,workflow,grounding"];
for (const [name, res] of Object.entries(allResults)) {
  if (!res.passK) continue;
  const fmt = v => (v === null || v === undefined) ? "" : v.toFixed(4);
  csvLines.push([
    name,
    MODEL_DISPLAY_NAMES[name] || name,
    fmt(res.passK.avg),
    fmt(res.passK.avgAttemptSuccessRate),
    fmt(res.quality?.avgWorkflow),
    fmt(res.quality?.avgGrounding)
  ].join(","));
}
fs.writeFileSync(outCsv, csvLines.join("\n"));

log(`\nResults saved → ${outJson}`);
log(`CSV → ${outCsv}`);
log(`\n✅ Done.`);
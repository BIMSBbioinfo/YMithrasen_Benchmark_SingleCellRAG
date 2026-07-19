/**
 * benchmark_runner.js
 *
 * Runs benchmark suites for the single-cell RAG evaluation.
 *
 * Suites:
 *   --suite internal   Pass@k, Workflow Success Rate, Grounding Rate (API models)
 *   --suite routing    Routing accuracy of ROLV domain classifier
 *   --suite qa         Normalized Gene-Turing Score
 *
 * Options:
 *   --k <number>       Number of attempts for Pass@k (default: 3)
 *   --model <names>    Comma-separated model keys to run, or "all" (default: all)
 *
 * Note: CLI-agent benchmarking (Gemini/Claude/Codex with/without recipe) now
 * lives in a separate benchmark project (coding_agent_benchmark), not here.
 */

import "dotenv/config";
import fs   from "fs";
import path from "path";
import { execFile, exec, execFileSync } from "child_process";
import { fileURLToPath }  from "url";
import { createReadStream } from "fs";
import readline from "readline";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ── Config ────────────────────────────────────────────────────────────────────

const CONFIG = {
  k:            parseInt(process.argv.find((a, i) => process.argv[i-1] === "--k") || "3"),
  suite:        process.argv.find((a, i) => process.argv[i-1] === "--suite") || "internal",
  model:        process.argv.find((a, i) => process.argv[i-1] === "--model") || "all",
  allowPartialPrompts: process.argv.includes("--allow-partial-prompts"),
  r_exec:       process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",
  dataset_dir:  (process.env.BENCH_DATASET_DIR ||
    "C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19"
  ).replace(/\\/g, "/"),
  results_dir:  path.join(__dirname, "results"),
  timeout_ms:   600000, // 10 min for complex pipelines with FindAllMarkers
  prompts_file: path.join(__dirname, "prompts/benchmark_prompts.jsonl"),
  routing_file: path.join(__dirname, "prompts/routing_queries.jsonl"),
  qa_file:      path.join(__dirname, "prompts/genomics_qa.jsonl"),
};

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 23);
fs.mkdirSync(CONFIG.results_dir, { recursive: true });

// ── Logging ───────────────────────────────────────────────────────────────────

function log(msg) {
  const line = `[${new Date().toISOString()}] ${msg}`;
  console.log(line);
}

// ── Load prompts ──────────────────────────────────────────────────────────────

function loadJSONL(file) {
  return fs.readFileSync(file, "utf8")
    .split("\n").filter(Boolean)
    .map(l => JSON.parse(l));
}

// ── Call baseline shim ────────────────────────────────────────────────────────

function callBaseline(model, prompt) {
  return new Promise((resolve) => {
    const payload = Buffer.from(JSON.stringify({ model, prompt })).toString("base64");
    const shimPath = path.join(__dirname, "scripts/baseline_shim.js");
    exec(
      `node "${shimPath}" "${payload}"`,
      { timeout: 60000 },
      (error, stdout) => {
        try {
          const parsed = JSON.parse(stdout?.trim().split("\n").filter(l => l.startsWith("{")).pop() || "{}");
          resolve(parsed);
        } catch {
          resolve({ result: null, error: error?.message || "parse error" });
        }
      }
    );
  });
}

// ── Call QA shim (raw text, no R extraction) ─────────────────────────────────

function callQA(model, prompt) {
  return new Promise((resolve) => {
    const shimPath = path.join(__dirname, "scripts/qa_shim.js");
    const payload  = Buffer.from(JSON.stringify({ model, prompt })).toString("base64");
    exec(
      `node "${shimPath}" "${payload}"`,
      { timeout: 60000 },
      (error, stdout) => {
        try {
          const line   = (stdout || "").trim().split("\n").filter(l => l.startsWith("{")).pop() || "{}";
          const parsed = JSON.parse(line);
          resolve({ result: parsed.result || null, error: parsed.error || error?.message || null });
        } catch {
          resolve({ result: null, error: error?.message || "parse error" });
        }
      }
    );
  });
}

// ── Call ROLV orchestrator shim ───────────────────────────────────────────────

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

// ── Task-specific validation ──────────────────────────────────────────────────
// Pass@k should mean "the task was actually accomplished," not merely
// "R didn't crash." This builds assertion code from a prompt's `validation`
// field and appends it after the generated code, before BENCH_EXEC_OK.
// A prompt with no `validation` field falls back to execution-only success
// (documented, not silently assumed) — see the `validated` flag in results.
function buildValidationR(validation = {}) {
  const checks = [];
  const explicitObj = validation.primary_object;
  // Only these fields reference the primary object — discovery is only needed
  // (and only counts as "this prompt has validation") when one of them is set.
  const needsPrimaryObject =
    validation.min_cells != null ||
    validation.min_features != null ||
    validation.required_reductions?.length ||
    validation.require_clusters ||
    validation.required_columns?.length;

  // Auto-discovery prelude: if the prompt doesn't require a specific variable
  // name, find whatever Seurat object exists in the global environment rather
  // than hard-requiring "pbmc". This measures task success, not naming
  // convention conformity. `.bench_obj` is used by all generic checks below;
  // if the prompt DOES specify primary_object, that name is used directly
  // instead (for prompts that legitimately require a specific output name,
  // e.g. "markers").
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

// ── Execute R code ────────────────────────────────────────────────────────────

const INSTALL_PREFIX = `
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
`.trim();

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
    INSTALL_PREFIX,
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
        // If no validation was defined for this prompt, success = execution only
        // (documented via `validated: false` in the result, not silently assumed).
        const validationSuccess = hasValidation ? (stdout || "").includes("BENCH_VALIDATION_OK") : true;
        const success = executionSuccess && validationSuccess;

        if (!success) {
          // Show the TAIL, not the head — R error messages print at the end of
          // stderr, after Seurat's verbose progress bars (log-normalization,
          // gene variance calc, etc). Slicing from the front was hiding the
          // actual error and only showing benign progress noise.
          log(`  [R exit] error=${error?.code} duration=${duration_ms}ms`);
          log(`  [R stderr] ${tail(stderr, 600)}`);
          log(`  [R stdout] ${tail(stdout, 400)}`);
        } else {
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        // Store full, untruncated stderr in results JSON — costs nothing and
        // preserves complete data for post-hoc debugging. Only the console
        // log preview above is shortened.
        resolve({
          success,
          executionSuccess,
          validationSuccess,
          validated: hasValidation,
          stderr: stderr || "",
          duration_ms,
        });
      }
    );
  });
}

function tail(text, n) {
  if (!text) return "(empty)";
  return text.length > n ? "…" + text.slice(-n) : text;
}

// ── Metrics ───────────────────────────────────────────────────────────────────
//
// Grounding rate is measured against a REAL function registry queried from R
// itself (installed package exports), not a hand-maintained guess-list. This
// means "grounded" genuinely means "resolves to an installed, real function or
// one the model defined itself" — not "happened to be on someone's list."

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

  // Write to a .R file and run by path — same pattern as executeRCode. Passing
  // a large multi-line, quote-heavy string via -e on Windows is a known source
  // of CreateProcess argument mangling; a file avoids that entirely.
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

// Detects functions the model defines itself (e.g. `run_qc <- function(x) {...}`)
// so legitimate local helpers aren't counted as hallucinated calls.
function extractDefinedFunctions(code) {
  const defined = new Set();
  const pattern = /\b([A-Za-z][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(/g;
  for (const m of code.matchAll(pattern)) defined.add(m[1]);
  return defined;
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
  // Strips # to end-of-line. Does not attempt to handle # inside string
  // literals (rare in generated Seurat code, and a false-negative here just
  // means a step is undercounted, not overcounted — safe failure direction).
  return code.split("\n").map(line => line.replace(/#.*$/, "")).join("\n");
}

// Workflow Success Rate is now driven by each prompt's own `expected_workflow`
// list (which steps THIS prompt actually requires), not by difficulty tier.
// Score = 70% step presence + 30% correct ordering among present steps.
// A prompt with no expected_workflow declared returns null (undeclared),
// not 0 — this is different from "failed to complete the workflow."
function calcWorkflowScore(code, expectedWorkflow = []) {
  if (!expectedWorkflow.length) return { score: null, presenceScore: null, orderScore: null, steps: [] };
  if (!code) return { score: 0, presenceScore: 0, orderScore: 0, steps: [] };

  code = stripRComments(code);

  const evaluated = expectedWorkflow.map(stepName => {
    const step = STEP_MAP.get(stepName);
    if (!step) return { name: stepName, present: false, position: -1 };
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

function wordOverlap(pred, ref) {
  const predWords = new Set(pred.toLowerCase().split(/\W+/).filter(Boolean));
  const refWords  = ref.toLowerCase().split(/\W+/).filter(Boolean);
  if (!refWords.length) return 0;
  return refWords.filter(w => predWords.has(w)).length / refWords.length;
}

// Detects whether k attempts actually differ (sampling diversity) vs producing
// near-identical code (which would make Pass@k ≈ Pass@1 regardless of k).
function normalizeCodeForComparison(code) {
  return (code || "").replace(/#.*$/gm, "").replace(/\s+/g, " ").trim();
}

function countUniqueGenerations(attempts) {
  const unique = new Set(attempts.map(a => normalizeCodeForComparison(a.code)).filter(Boolean));
  return unique.size;
}

// ── Pass@k suite ──────────────────────────────────────────────────────────────

async function runPassK(modelName, callerFn, prompts, k) {
  log(`\n[Pass@k] model=${modelName}  k=${k}  prompts=${prompts.length}`);
  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < k; i++) {
      const res  = await callerFn(p.prompt);
      // baseline_shim.js already extracts R code from fences before returning;
      // re-running extraction heuristics on that clean code was rejecting valid
      // short snippets (e.g. a single VlnPlot() call) that don't match the
      // library(/<-/Read10X fallback markers extractRFromText looks for.
      const code = res.alreadyExtracted
        ? (res.result || null)
        : extractRFromText(res.result || "");
      if (!code) {
        attempts.push({ passed: false, error: "no_r_code_extracted", model_error: res.error || null });
        continue;
      }
      const exec = await executeRCode(code, CONFIG.dataset_dir, p.validation || {});
      attempts.push({
        passed:            exec.success,
        code,
        executionSuccess:  exec.executionSuccess,
        validationSuccess: exec.validationSuccess,
        validated:         exec.validated,
        error:             exec.success ? null : (exec.stderr || "unknown error").trim(),
        duration_ms:       exec.duration_ms,
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
  // If no code fence found, treat the whole text as R code
  // (rolv-raw returns plain R code without fences)
  if (text.includes("library(") || text.includes("<-") || text.includes("Read10X")) {
    return text;
  }
  return null;
}

// ── Code quality suite ────────────────────────────────────────────────────────
// Reuses the code already generated during Pass@k rather than calling the model
// again — no extra API cost. Scores EVERY attempt (not just the first), since
// Pass@k considers all k attempts and workflow/grounding should describe the
// same population of generations, not a single arbitrary one.
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

    // null-safe averaging: workflow.score is null when a prompt declares no
    // expected_workflow, and grounding.rate is null when the registry failed to build
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

// ── Routing suite ─────────────────────────────────────────────────────────────

async function runRouting() {
  log(`\n[Routing] evaluating ROLV classifier`);
  const queries = loadJSONL(CONFIG.routing_file);
  const { default: axios } = await import("axios").catch(() => ({ default: null }));

  // Call ROLV classifier endpoint directly
  const correct = { code:0, pubmed:0, geo:0, vision:0, general:0 };
  const total   = { code:0, pubmed:0, geo:0, vision:0, general:0 };
  const tp = {}, fp = {}, fn_ = {};
  const domains = ["code","pubmed","geo","vision","general"];
  domains.forEach(d => { tp[d]=0; fp[d]=0; fn_[d]=0; });

  for (const q of queries) {
    const res = await callRolv(q.query);
    // Extract domain from ROLV response metadata if available
    // For now, score based on whether we got a result at all
    total[q.label] = (total[q.label] || 0) + 1;
  }

  log(`  routing evaluation requires ROLV classifier endpoint — see scripts/routing_eval.js`);
  return {};
}

// ── QA suite ──────────────────────────────────────────────────────────────────

async function runQA(modelName, callerFn, qaPairs) {
  log(`\n[QA / Gene-Turing] model=${modelName}  questions=${qaPairs.length}`);
  const scores = [];

  // All QA calls go through qa_shim.js which returns raw text — no R code extraction.
  // baseline_shim.js strips non-R content so can't be used for factual QA.
  // rolv/rolv-nopreexec use "claude" model key; groq uses "groq" model key.
  const qaModel = modelName === "groq" ? "groq" : "claude";

  for (const qa of qaPairs) {
    const res = await callQA(qaModel, qa.question);

    // Extract text — handle both { result: "..." } and plain string responses
    let answer = "";
    if (typeof res === "string") answer = res;
    else if (res && typeof res.result === "string") answer = res.result;
    else if (res && typeof res === "object") answer = JSON.stringify(res);

    const score = wordOverlap(answer, qa.answer);
    scores.push({ id: qa.id, question: qa.question, expected: qa.answer, got: answer.slice(0, 150), score });
    log(`  [${qa.id}] score=${score.toFixed(2)}  expected="${qa.answer}"  got="${answer.slice(0, 80)}"`);
  }

  const avg = scores.reduce((s, r) => s + r.score, 0) / scores.length;
  log(`  avg Gene-Turing Score = ${avg.toFixed(3)}`);
  return { avg, scores };
}

// ── Main ──────────────────────────────────────────────────────────────────────

// ── Preflight checks ──────────────────────────────────────────────────────────

// Maps each model key to the env vars it actually needs. Preflight only
// requires/smoke-tests keys for models actually selected via --model, so a
// downed ROLV app (or a missing OpenRouter key, etc.) never blocks a run that
// doesn't touch that provider.
const MODEL_REQUIRED_KEYS = {
  "rolv":                ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "rolv-nopreexec":      ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "rolv-raw":            ["ROLV_API_URL", "ROLV_JWT_TOKEN"],
  "groq":                ["GROQ_API_KEY"],
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
  GROQ_API_KEY:       "Groq key",
  ANTHROPIC_API_KEY:  "Anthropic key",
};

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

  // Catch typos (e.g. "claude-sonet-5") before they silently drop a model
  // from the comparison — a missing model in a final thesis run is worse
  // than an aborted run.
  const unknownModels = selectedModels.filter(
    m => !Object.prototype.hasOwnProperty.call(MODEL_REQUIRED_KEYS, m)
  );
  if (unknownModels.length) {
    errors.push(
      `Unknown model key(s): ${unknownModels.join(", ")} — check spelling against: ${Object.keys(MODEL_REQUIRED_KEYS).join(", ")}`
    );
  }

  // 1. Required env vars — union of keys needed by the selected (recognized) models only
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

  // 2. JWT expiry check — only if a ROLV model is in scope
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

  // 3. R executable exists
  if (!fs.existsSync(CONFIG.r_exec)) {
    errors.push(`R executable not found: ${CONFIG.r_exec}`);
  } else {
    log(`  ✅ R executable found: ${CONFIG.r_exec}`);
  }

  // 4. Dataset directory exists
  if (!fs.existsSync(CONFIG.dataset_dir)) {
    errors.push(`Dataset directory not found: ${CONFIG.dataset_dir}`);
  } else {
    log(`  ✅ Dataset directory found: ${CONFIG.dataset_dir}`);
  }

  // 5. Prompts file exists, and (for internal suite) check validation coverage.
  // Missing validation/expected_workflow is a HARD ERROR by default — the final
  // Pass@3 run should not silently mix task-validated and execution-only results.
  // Pass --allow-partial-prompts for quick exploratory/smoke-test runs where
  // not every prompt is annotated yet.
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
      } catch (e) {
        warnings.push(`Could not pre-check prompt validation coverage: ${e.message}`);
      }
    }
  }

  // 6. Live API smoke tests — runs even if static checks passed, catches bad/expired keys.
  // Only smoke-tests models actually in scope for this run.
  if (!errors.length) {
    log("\n  Running live API smoke tests (1 call each, selected models only)...");
    const SMOKE_PROMPT = "Write a single line of R code that assigns the number 1 to a variable called x.";

    const ALL_SMOKE_TESTS = {
      "gpt-5.5":            () => callBaseline("gpt-5.5", SMOKE_PROMPT),
      "gpt-5.4-mini":       () => callBaseline("gpt-5.4-mini", SMOKE_PROMPT),
      "gemini-openrouter":  () => callBaseline("gemini-openrouter", SMOKE_PROMPT),
      "gemini-flash-lite":  () => callBaseline("gemini-flash-lite", SMOKE_PROMPT),
      "groq":               () => callBaseline("groq", SMOKE_PROMPT),
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

// Build the R function registry once — used by every model's grounding score.
// This is a real query against installed packages, not a hardcoded guess-list.
log("\nBuilding R function registry for grounding rate...");
const R_FUNCTION_REGISTRY = getRFunctionRegistry();

const prompts = loadJSONL(CONFIG.prompts_file);
const qaPairs = CONFIG.suite === "qa" ? loadJSONL(CONFIG.qa_file) : [];

const API_MODELS = {
  rolv:                (p) => callRolv(p),
  "rolv-nopreexec":    (p) => callRolv(p, "orchestrator_shim_nopreexec.js"),
  "rolv-raw":          (p) => callBaseline("rolv-raw", p),
  groq:                (p) => callBaseline("groq", p),
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
  "groq":              "GPT-OSS 20B via Groq",
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

if (CONFIG.suite === "qa" || CONFIG.suite === "all") {
  log(`\n══ QA Suite ══`);
  for (const [name, fn] of Object.entries(API_MODELS)) {
    const qa = await runQA(name, fn, qaPairs);
    allResults[name] = { ...(allResults[name] || {}), qa };
  }
}

if (CONFIG.suite === "routing") {
  await runRouting();
}

// ── Summary ───────────────────────────────────────────────────────────────────

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
  "groq":                { temperature: 0.6 },
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
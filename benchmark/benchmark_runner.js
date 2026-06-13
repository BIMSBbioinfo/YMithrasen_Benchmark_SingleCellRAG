/**
 * benchmark_runner.js
 *
 * Runs benchmark suites for the single-cell RAG evaluation.
 *
 * Suites:
 *   --suite internal   Pass@k, Workflow Success Rate, Grounding Rate (API models)
 *   --suite routing    Routing accuracy of ROLV domain classifier
 *   --suite qa         Normalized Gene-Turing Score
 *   --suite cli        CLI agents with and without recipe (Gemini, Claude, Codex)
 *
 * Options:
 *   --k <number>       Number of attempts for Pass@k (default: 3)
 *   --agent <name>     For CLI suite: gemini | claude | codex | all (default: all)
 */

import "dotenv/config";
import fs   from "fs";
import path from "path";
import { execFile, exec } from "child_process";
import { fileURLToPath }  from "url";
import { createReadStream } from "fs";
import readline from "readline";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ── Config ────────────────────────────────────────────────────────────────────

const CONFIG = {
  k:            parseInt(process.argv.find((a, i) => process.argv[i-1] === "--k") || "3"),
  suite:        process.argv.find((a, i) => process.argv[i-1] === "--suite") || "internal",
  agent:        process.argv.find((a, i) => process.argv[i-1] === "--agent") || "all",
  r_exec:       process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",
  dataset_dir:  (process.env.BENCH_DATASET_DIR ||
    "C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19"
  ).replace(/\\/g, "/"),
  recipe_path:  path.join(__dirname, "../recipes/scrna-seurat-pbmc/recipe.md").replace(/\\/g, "/"),
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

// ── Call CLI runner ───────────────────────────────────────────────────────────

function callCLI(agent, prompt, recipeEnabled, workDir) {
  return new Promise((resolve) => {
    const recipe_path = recipeEnabled ? CONFIG.recipe_path : null;
    const payload = Buffer.from(JSON.stringify({
      agent, prompt, recipe_path, work_dir: workDir
    })).toString("base64");
    const runnerPath = path.join(__dirname, "scripts/cli_runner.js");
    exec(
      `node "${runnerPath}" "${payload}"`,
      { timeout: 120000 },
      (error, stdout) => {
        try {
          const parsed = JSON.parse(stdout?.trim() || "{}");
          resolve(parsed);
        } catch {
          resolve({ result: null, error: error?.message || "parse error" });
        }
      }
    );
  });
}

// ── Execute R code ────────────────────────────────────────────────────────────

const INSTALL_PREFIX = `
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
`.trim();

function executeRCode(rCode, datasetDir) {
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

  const fullScript = [
    `setwd("${datasetDir.replace(/\\/g, "/")}")`,
    INSTALL_PREFIX,
    `tryCatch({`,
    `  suppressMessages(suppressWarnings({`,
    rCode,
    `  }))`,
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
        const duration_ms = Date.now() - start;
        const success = !error && stdout.includes("BENCH_EXEC_OK");
        if (!success) {
          log(`  [R exit] error=${error?.code} duration=${duration_ms}ms`);
          log(`  [R stderr] ${(stderr || "(empty)").slice(0, 400)}`);
          log(`  [R stdout] ${(stdout || "(empty)").slice(0, 400)}`);
        } else {
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        resolve({ success, stderr: (stderr || "").slice(0, 400), duration_ms });
      }
    );
  });
}

// ── Metrics ───────────────────────────────────────────────────────────────────

const SEURAT_FUNCTIONS = new Set([
  "Read10X","ReadMtx","CreateSeuratObject","PercentageFeatureSet","VlnPlot",
  "FeatureScatter","subset","NormalizeData","FindVariableFeatures","ScaleData",
  "RunPCA","ElbowPlot","DimHeatmap","VizDimLoadings","FindNeighbors",
  "FindClusters","RunUMAP","RunTSNE","DimPlot","FeaturePlot","DoHeatmap",
  "FindAllMarkers","FindMarkers","RenameIdents","WhichCells","SCTransform",
  "IntegrationAnchorSet","FindIntegrationAnchors","IntegrateData",
  "IntegrateLayers","CellCycleScoring","AddMetaData","Embeddings",
  "GetAssayData","LayerData","JoinLayers","SplitObject","merge",
  "paramSweep","summarizeSweep","find.pK","modelHomotypic","doubletFinder",
  "library","require","print","head","cat","paste","paste0",
  "ggplot","aes","geom_point","geom_line","xlab","ylab","theme",
  "NoLegend","LabelPoints","VariableFeaturePlot",
  "filter","select","mutate","group_by","arrange","slice_head","ungroup",
  "rownames","colnames","names","levels","nrow","ncol","round",
  "saveRDS","readRDS","setwd","getwd"
]);

const PIPELINE_STEPS = [
  { name: "data_loading",    patterns: ["Read10X","ReadMtx","read.table","read.csv"] },
  { name: "seurat_object",   patterns: ["CreateSeuratObject"] },
  { name: "qc",              patterns: ["PercentageFeatureSet","VlnPlot","FeatureScatter"] },
  { name: "normalisation",   patterns: ["NormalizeData","SCTransform"] },
  { name: "variable_feats",  patterns: ["FindVariableFeatures"] },
  { name: "scaling",         patterns: ["ScaleData"] },
  { name: "pca",             patterns: ["RunPCA"] },
  { name: "graph",           patterns: ["FindNeighbors"] },
  { name: "clustering",      patterns: ["FindClusters"] },
  { name: "embedding",       patterns: ["RunUMAP","RunTSNE"] },
];

function calcWorkflowScore(code) {
  if (!code) return 0;
  const found = PIPELINE_STEPS.filter(s => s.patterns.some(p => code.includes(p)));
  return found.length / PIPELINE_STEPS.length;
}

function calcGroundingRate(code) {
  if (!code) return 0;
  const calls = [...code.matchAll(/\b([A-Za-z][A-Za-z0-9_.]*)\s*\(/g)].map(m => m[1]);
  if (!calls.length) return 1;
  const grounded = calls.filter(c => SEURAT_FUNCTIONS.has(c));
  return grounded.length / calls.length;
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

// ── Pass@k suite ──────────────────────────────────────────────────────────────

async function runPassK(modelName, callerFn, prompts, k) {
  log(`\n[Pass@k] model=${modelName}  k=${k}  prompts=${prompts.length}`);
  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < k; i++) {
      const res  = await callerFn(p.prompt);
      const code = extractRFromText(res.result || "");
      if (!code) { attempts.push({ passed: false }); continue; }
      const exec = await executeRCode(code, CONFIG.dataset_dir);
      attempts.push({ passed: exec.success, code });
    }
    const c = attempts.filter(a => a.passed).length;
    const score = passAtK(k, c, k);
    log(`  [${p.id}] pass@${k}=${score.toFixed(3)}  (${c}/${k} passed)`);
    results.push({ id: p.id, score, attempts });
  }

  const avg = results.reduce((s, r) => s + r.score, 0) / results.length;
  log(`  avg Pass@${k} = ${avg.toFixed(3)}`);
  return { avg, results };
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

async function runCodeQuality(modelName, callerFn, prompts) {
  log(`\n[CodeQuality] model=${modelName}`);
  const scores = [];

  for (const p of prompts) {
    const res  = await callerFn(p.prompt);
    const code = res.result || extractRFromText(res.result || "");
    scores.push({
      id:        p.id,
      workflow:  calcWorkflowScore(code),
      grounding: calcGroundingRate(code)
    });
  }

  const avgWorkflow  = scores.reduce((s, r) => s + r.workflow,  0) / scores.length;
  const avgGrounding = scores.reduce((s, r) => s + r.grounding, 0) / scores.length;
  log(`  avg workflow success rate: ${avgWorkflow.toFixed(3)}`);
  log(`  avg grounding rate:        ${avgGrounding.toFixed(3)}`);
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

// ── CLI suite ─────────────────────────────────────────────────────────────────

async function runCLISuite(agents, prompts, k) {
  log(`\n══ CLI Suite ══`);
  const results = {};

  for (const agent of agents) {
    for (const withRecipe of [true, false]) {
      const condition = `${agent}_${withRecipe ? "with_recipe" : "vanilla"}`;
      log(`\n[CLI] agent=${agent}  recipe=${withRecipe}  k=${k}  prompts=${prompts.length}`);

      const tmpDir = path.join(CONFIG.dataset_dir, "..", `cli_tmp_${agent}`);
      fs.mkdirSync(tmpDir, { recursive: true });

      const callerFn = (prompt) => callCLI(agent, prompt, withRecipe, tmpDir);
      const passK    = await runPassK(condition, callerFn, prompts, k);
      const quality  = await runCodeQuality(condition, callerFn, prompts);

      results[condition] = { passK, quality };
    }
  }

  return results;
}

// ── Main ──────────────────────────────────────────────────────────────────────

log("═".repeat(60));
log(`ROLV Benchmark  run_id=${run_id}`);
log(`k=${CONFIG.k}  suite=${CONFIG.suite}`);
log("═".repeat(60));

const prompts = loadJSONL(CONFIG.prompts_file);
const qaPairs = CONFIG.suite === "qa" ? loadJSONL(CONFIG.qa_file) : [];

const API_MODELS = {
  rolv:           (p) => callRolv(p),                          // full RAG + pre-exec planner
  "rolv-nopreexec": (p) => callRolv(p, "orchestrator_shim_nopreexec.js"), // RAG only, no pre-exec
  "rolv-raw": (p) => callBaseline("rolv-raw", p), // direct API — vanilla Claude via @out-of-domain
  groq:       (p) => callBaseline("groq", p),
  // gemini: excluded — free tier quota exhausted
};

const allResults = {};

if (CONFIG.suite === "internal" || CONFIG.suite === "all") {
  log(`\n══ Internal Suite ══`);
  for (const [name, fn] of Object.entries(API_MODELS)) {
    log(`\n── ${name} ──`);
    const passK   = await runPassK(name, fn, prompts, CONFIG.k);
    const quality = await runCodeQuality(name, fn, prompts);
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

if (CONFIG.suite === "cli" || CONFIG.suite === "all") {
  const agentArg = CONFIG.agent === "all"
    ? ["gemini", "claude", "codex"]
    : [CONFIG.agent];
  const cliResults = await runCLISuite(agentArg, prompts, CONFIG.k);
  Object.assign(allResults, cliResults);
}

// ── Summary ───────────────────────────────────────────────────────────────────

log("\n" + "─".repeat(68));
log("BENCHMARK SUMMARY");
log("─".repeat(68));

if (allResults && Object.keys(allResults).length) {
  log("\nINTERNAL METRICS");
  log(`Model              Pass@k     Workflow    Grounding`);
  log("─".repeat(50));
  for (const [name, res] of Object.entries(allResults)) {
    if (!res.passK) continue;
    const p = (res.passK.avg * 100).toFixed(1) + "%";
    const w = (res.quality?.avgWorkflow  * 100 || 0).toFixed(1) + "%";
    const g = (res.quality?.avgGrounding * 100 || 0).toFixed(1) + "%";
    log(`${name.padEnd(18)} ${p.padEnd(10)} ${w.padEnd(12)} ${g}`);
  }
}

// Save results
const outJson = path.join(CONFIG.results_dir, `results_${run_id}.json`);
const outCsv  = path.join(CONFIG.results_dir, `results_${run_id}_summary.csv`);
fs.writeFileSync(outJson, JSON.stringify({ run_id, config: CONFIG, results: allResults }, null, 2));

// CSV
const csvLines = ["model,pass_k,workflow,grounding"];
for (const [name, res] of Object.entries(allResults)) {
  if (!res.passK) continue;
  csvLines.push([
    name,
    res.passK.avg.toFixed(4),
    (res.quality?.avgWorkflow  || 0).toFixed(4),
    (res.quality?.avgGrounding || 0).toFixed(4)
  ].join(","));
}
fs.writeFileSync(outCsv, csvLines.join("\n"));

log(`\nResults saved → ${outJson}`);
log(`CSV → ${outCsv}`);
log(`\n✅ Done.`);
/**
 * cli_benchmark.js
 *
 * Benchmarks Gemini CLI with and without the scRNA-seq recipe.
 *
 * With recipe:    writes recipe.md as GEMINI.md in the work dir so Gemini CLI
 *                 picks it up automatically as project context.
 * Without recipe: no GEMINI.md — agent works from scratch.
 *
 * Usage:
 *   node cli_benchmark.js --k 1
 *   node cli_benchmark.js --k 3
 */

import "dotenv/config";
import fs        from "fs";
import path      from "path";
import os        from "os";
import { exec, execFile } from "child_process";
import { fileURLToPath }  from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// ── Config ────────────────────────────────────────────────────────────────────

const K           = parseInt(process.argv.find((a,i) => process.argv[i-1]==="--k") || "1");
const R_EXEC      = (process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe");
const DATASET_DIR = (process.env.BENCH_DATASET_DIR ||
  "C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19"
).replace(/\\/g, "/");
const RECIPE_PATH = path.join(__dirname, "../recipes/scrna-seurat-pbmc/recipe.md");
const PROMPTS_FILE = path.join(__dirname, "prompts/benchmark_prompts.jsonl");
const RESULTS_DIR  = path.join(__dirname, "results");
const TIMEOUT_MS   = 180000; // 3 min per CLI call

fs.mkdirSync(RESULTS_DIR, { recursive: true });

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 23);

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

// ── Extract R code ────────────────────────────────────────────────────────────

function extractRCode(text) {
  const match = text.match(/```(?:\{r\}|r)\n?([\s\S]*?)```/);
  return match ? match[1].trim() : null;
}

// ── Call Gemini CLI ───────────────────────────────────────────────────────────

function callGeminiCLI(prompt, workDir) {
  return new Promise((resolve) => {
    const systemContext = `You are a bioinformatics assistant. When given a task involving R and single-cell RNA-seq analysis, respond with ONLY executable R code.
Use R Markdown code fence format: \`\`\`{r}
...code...
\`\`\`
Do NOT include ANY explanatory text outside the code block.
The data files are in the current working directory. Use Read10X(data.dir = ".") to load 10x data.`;

    const fullPrompt = `${systemContext}\n\nTask:\n${prompt}`;

    exec(
      `gemini -p "${fullPrompt.replace(/"/g, '\\"').replace(/\n/g, ' ')}" --output-format text --yolo`,
      { timeout: TIMEOUT_MS, cwd: workDir },
      (error, stdout, stderr) => {
        if (error) {
          log(`  [gemini-cli] error: ${error.message.slice(0, 200)}`);
          resolve({ result: null, error: error.message });
          return;
        }
        const code = extractRCode(stdout);
        resolve({ result: code, raw: stdout, error: code ? null : "No R code block found" });
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

function executeRCode(rCode) {
  const tmpId   = Math.random().toString(36).slice(2, 8);
  const tmpDir  = path.join(RESULTS_DIR, `cli_rexec_${tmpId}`);
  fs.mkdirSync(tmpDir, { recursive: true });
  const scriptPath = path.join(tmpDir, "bench.R");

  // Normalise paths
  rCode = rCode
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']\.\/hg19\/?["']\s*\)/gi, 'Read10X(data.dir = ".")')
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']hg19\/?["']\s*\)/gi, 'Read10X(data.dir = ".")')
    .replace(/["']path\/to\/[^"']*["']/gi, '"."');

  const fullScript = [
    `setwd("${DATASET_DIR}")`,
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

  return new Promise((resolve) => {
    const start = Date.now();
    execFile(R_EXEC, ["--vanilla", scriptPath], { timeout: 120000 },
      (error, stdout, stderr) => {
        const duration_ms = Date.now() - start;
        const success = !error && stdout.includes("BENCH_EXEC_OK");
        if (!success) {
          log(`  [R stderr] ${(stderr || "(empty)").slice(0, 300)}`);
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
  "RunPCA","ElbowPlot","DimHeatmap","FindNeighbors","FindClusters","RunUMAP",
  "RunTSNE","DimPlot","FeaturePlot","DoHeatmap","FindAllMarkers","FindMarkers",
  "RenameIdents","SCTransform","JoinLayers","paramSweep","summarizeSweep",
  "find.pK","modelHomotypic","doubletFinder","library","require",
  "ggplot","aes","geom_point","geom_line","xlab","ylab","theme","NoLegend",
  "filter","select","mutate","group_by","arrange","slice_head","ungroup",
  "rownames","colnames","names","levels","nrow","ncol","round","saveRDS","readRDS"
]);

const PIPELINE_STEPS = [
  { name: "data_loading",   patterns: ["Read10X","ReadMtx","read.table","read.csv"] },
  { name: "seurat_object",  patterns: ["CreateSeuratObject"] },
  { name: "qc",             patterns: ["PercentageFeatureSet","VlnPlot","FeatureScatter"] },
  { name: "normalisation",  patterns: ["NormalizeData","SCTransform"] },
  { name: "variable_feats", patterns: ["FindVariableFeatures"] },
  { name: "scaling",        patterns: ["ScaleData"] },
  { name: "pca",            patterns: ["RunPCA"] },
  { name: "graph",          patterns: ["FindNeighbors"] },
  { name: "clustering",     patterns: ["FindClusters"] },
  { name: "embedding",      patterns: ["RunUMAP","RunTSNE"] },
];

function calcWorkflowScore(code) {
  if (!code) return 0;
  return PIPELINE_STEPS.filter(s => s.patterns.some(p => code.includes(p))).length / PIPELINE_STEPS.length;
}

function calcGroundingRate(code) {
  if (!code) return 0;
  const calls = [...code.matchAll(/\b([A-Za-z][A-Za-z0-9_.]*)\s*\(/g)].map(m => m[1]);
  if (!calls.length) return 1;
  return calls.filter(c => SEURAT_FUNCTIONS.has(c)).length / calls.length;
}

function passAtK(n, c, k) {
  if (n - c < k) return 1.0;
  let num = 1, den = 1;
  for (let i = 0; i < k; i++) { num *= (n - c - i); den *= (n - i); }
  return 1 - num / den;
}

// ── Run one condition ─────────────────────────────────────────────────────────

async function runCondition(label, prompts, withRecipe) {
  log(`\n${"─".repeat(60)}`);
  log(`Condition: ${label}  (recipe=${withRecipe}  k=${K})`);
  log("─".repeat(60));

  // Create a shared work dir for this condition
  const workDir = path.join(RESULTS_DIR, `cli_workdir_${label.replace(/\s/g,"_")}`);
  fs.mkdirSync(workDir, { recursive: true });

  // Write or remove GEMINI.md
  const geminiMdPath = path.join(workDir, "GEMINI.md");
  if (withRecipe) {
    const recipeContent = fs.readFileSync(RECIPE_PATH, "utf8");
    fs.writeFileSync(geminiMdPath, `# scRNA-seq Workflow Recipe\n\n${recipeContent}`);
    log(`  [recipe] GEMINI.md written (${recipeContent.length} chars)`);
  } else {
    if (fs.existsSync(geminiMdPath)) fs.unlinkSync(geminiMdPath);
    log(`  [recipe] no GEMINI.md — vanilla mode`);
  }

  const results = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < K; i++) {
      log(`  [${p.id}] attempt ${i+1}/${K} ...`);
      const res  = await callGeminiCLI(p.prompt, workDir);
      if (!res.result) {
        log(`  [${p.id}] no code extracted`);
        attempts.push({ passed: false, code: null });
        continue;
      }
      const exec = await executeRCode(res.result);
      log(`  [${p.id}] R execution: ${exec.success ? "PASS" : "FAIL"} (${exec.duration_ms}ms)`);
      attempts.push({ passed: exec.success, code: res.result });
    }

    const c     = attempts.filter(a => a.passed).length;
    const score = passAtK(K, c, K);
    const bestCode = attempts.find(a => a.code)?.code || null;

    results.push({
      id:        p.id,
      score,
      c, n: K,
      workflow:  calcWorkflowScore(bestCode),
      grounding: calcGroundingRate(bestCode),
      attempts
    });

    log(`  [${p.id}] pass@${K}=${score.toFixed(3)}  (${c}/${K} passed)`);
  }

  const avgPassK     = results.reduce((s,r) => s + r.score, 0) / results.length;
  const avgWorkflow  = results.reduce((s,r) => s + r.workflow, 0) / results.length;
  const avgGrounding = results.reduce((s,r) => s + r.grounding, 0) / results.length;

  log(`\n  avg Pass@${K}           = ${avgPassK.toFixed(3)}`);
  log(`  avg workflow score     = ${avgWorkflow.toFixed(3)}`);
  log(`  avg grounding rate     = ${avgGrounding.toFixed(3)}`);

  return { label, withRecipe, avgPassK, avgWorkflow, avgGrounding, results };
}

// ── Main ──────────────────────────────────────────────────────────────────────

log("═".repeat(68));
log(`CLI Benchmark  run_id=${run_id}  k=${K}`);
log("═".repeat(68));

const prompts = loadJSONL(PROMPTS_FILE);

const withRecipe    = await runCondition("gemini_with_recipe",    prompts, true);
const withoutRecipe = await runCondition("gemini_without_recipe", prompts, false);

// ── Summary ───────────────────────────────────────────────────────────────────

log("\n" + "═".repeat(68));
log("CLI BENCHMARK SUMMARY");
log("═".repeat(68));
log(`\nCondition                  Pass@${K}    Workflow    Grounding`);
log("─".repeat(56));

for (const c of [withRecipe, withoutRecipe]) {
  const p = (c.avgPassK     * 100).toFixed(1) + "%";
  const w = (c.avgWorkflow  * 100).toFixed(1) + "%";
  const g = (c.avgGrounding * 100).toFixed(1) + "%";
  log(`${c.label.padEnd(26)} ${p.padEnd(10)} ${w.padEnd(12)} ${g}`);
}

const deltaPass = ((withRecipe.avgPassK - withoutRecipe.avgPassK) * 100).toFixed(1);
const deltaWork = ((withRecipe.avgWorkflow - withoutRecipe.avgWorkflow) * 100).toFixed(1);
const deltaGround = ((withRecipe.avgGrounding - withoutRecipe.avgGrounding) * 100).toFixed(1);
log(`\nDelta (with - without):    ${deltaPass > 0 ? "+" : ""}${deltaPass}%     ${deltaWork > 0 ? "+" : ""}${deltaWork}%       ${deltaGround > 0 ? "+" : ""}${deltaGround}%`);

// Save results
const output = {
  run_id,
  k: K,
  conditions: [withRecipe, withoutRecipe],
  delta: { pass_k: deltaPass, workflow: deltaWork, grounding: deltaGround }
};

const outJson = path.join(RESULTS_DIR, `cli_results_${run_id}.json`);
const outCsv  = path.join(RESULTS_DIR, `cli_results_${run_id}.csv`);

fs.writeFileSync(outJson, JSON.stringify(output, null, 2));

const csvLines = ["condition,with_recipe,pass_k,workflow,grounding"];
for (const c of [withRecipe, withoutRecipe]) {
  csvLines.push([c.label, c.withRecipe, c.avgPassK.toFixed(4), c.avgWorkflow.toFixed(4), c.avgGrounding.toFixed(4)].join(","));
}
fs.writeFileSync(outCsv, csvLines.join("\n"));

log(`\nResults → ${outJson}`);
log(`CSV     → ${outCsv}`);
log(`\n✅ Done.`);

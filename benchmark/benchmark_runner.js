/**
 * benchmark_runner.js
 *
 * Main benchmark harness for ROLV Agentic RAG system.
 * Integrates directly with the existing orchestrator.js pipeline.
 *
 * Usage:
 *   node benchmark_runner.js                      # run all suites
 *   node benchmark_runner.js --suite internal     # only internal metrics
 *   node benchmark_runner.js --suite routing      # only routing accuracy
 *   node benchmark_runner.js --suite qa           # only Genomics QA
 *   node benchmark_runner.js --k 3                # set Pass@k value (default 3)
 *   BENCH_DATASET_DIR=/path/to/data node benchmark_runner.js
 */

import fs from "fs";
import path from "path";
import { exec, execFile } from "child_process";
import { promisify } from "util";
import { fileURLToPath } from "url";

const execAsync = promisify(exec);
const __dirname = path.dirname(fileURLToPath(import.meta.url));

// CONFIG
const ARGS = Object.fromEntries(
  process.argv.slice(2).reduce((acc, v, i, arr) => {
    if (v.startsWith("--")) acc.push([v.slice(2), arr[i + 1]]);
    return acc;
  }, [])
);

const CONFIG = {
  k:            parseInt(ARGS.k || process.env.BENCH_K || "3"),
  timeout_ms:   120_000,
  results_dir:  path.join(__dirname, "results"),
  prompts_file: path.join(__dirname, "prompts", "benchmark_prompts.jsonl"),
  routing_file: path.join(__dirname, "prompts", "routing_queries.jsonl"),
  qa_file:      path.join(__dirname, "prompts", "genomics_qa.jsonl"),
  r_exec:       process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",
  dataset_dir:  (process.env.BENCH_DATASET_DIR || "C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19").replace(/\\/g, "/"),
  run_id:       new Date().toISOString().replace(/[:.]/g, "-"),
};

// LOGGING
fs.mkdirSync(CONFIG.results_dir, { recursive: true });
const logPath = path.join(CONFIG.results_dir, `run_${CONFIG.run_id}.log`);

function log(msg) {
  const line = `[${new Date().toISOString()}] ${msg}`;
  console.log(line);
  fs.appendFileSync(logPath, line + "\n");
}

// RESULT STORE

const results = { run_id: CONFIG.run_id, internal: {}, routing: null, qa: {}, summary: {} };

function saveResults() {
  const outPath = path.join(CONFIG.results_dir, `results_${CONFIG.run_id}.json`);
  fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
  log(`Results saved → ${outPath}`);
  return outPath;
}

// MODEL ADAPTERS

/**
 * The shim imports orchestrator.js,
 * calls orchestrator(prompt), and writes { result } to stdout as JSON.
 * This allows us to reuse the exact same orchestrator logic + prompt engineering as the main app,
 * while still running the benchmark as a standalone script.
 */
async function callOrchestrator(prompt, history = []) {
  const shimPath = path.join(__dirname, "scripts", "orchestrator_shim.js");
  if (!fs.existsSync(shimPath)) {
    log(`[orchestrator] shim missing at ${shimPath} — see README`);
    return null;
  }
  const payload = Buffer.from(JSON.stringify({ prompt, history })).toString("base64");
  try {
    const { stdout } = await execAsync(
      `node "${shimPath}" "${payload}"`,
      { timeout: CONFIG.timeout_ms }
    );
    // Find the last line that looks like JSON — [orch-log] lines come first,
    // the actual {"result":...} line is always last
    const jsonLine = stdout.split("\n").map(l => l.trim()).filter(l => l.startsWith("{")).pop();
    if (!jsonLine) { log(`[orchestrator] no JSON line in output`); return null; }
    return JSON.parse(jsonLine).result || null;
  } catch (err) {
    log(`[orchestrator] error: ${err.message.slice(0, 200)}`);
    return null;
  }
}

async function callBaseline(model, prompt) {
  const shimPath = path.join(__dirname, "scripts", "baseline_shim.js");
  if (!fs.existsSync(shimPath)) {
    log(`[baseline:${model}] shim missing — skipping`);
    return null;
  }
  const payload = Buffer.from(JSON.stringify({ model, prompt })).toString("base64");
  try {
    const { stdout } = await execAsync(
      `node "${shimPath}" "${payload}"`,
      { timeout: CONFIG.timeout_ms }
    );
    return JSON.parse(stdout.trim()).result || null;
  } catch (err) {
    log(`[baseline:${model}] error: ${err.message.slice(0, 200)}`);
    return null;
  }
}

const MODELS = {
  rolv:   (p) => callOrchestrator(p),
  gpt4o:  (p) => callBaseline("gpt-4o", p),
  gemini: (p) => callBaseline("gemini-2.0-flash", p),
};

// R CODE EXTRACTOR
// Handles all fence variants the app uses (from mainChat.js system prompt):
//   ```{r}          — standard R Markdown
//   ```{r, echo=FALSE}  — R Markdown with options
//   ```r            — plain r fence (some models output this)
function extractRCode(text) {
  if (!text) return null;
  // Priority order matches what mainChat.js enforces via its system prompt:
  // "Always return R code inside executable R Markdown chunks: use ```{r}```"
  const patterns = [
    /```\{r[^}]*\}\s*\n([\s\S]*?)```/i,   // ```{r}  or  ```{r, echo=FALSE}
    /```r\s*\n([\s\S]*?)```/i,             // ```r (plain fence, some baselines output this)
    /```R\s*\n([\s\S]*?)```/,              // ```R (uppercase, seen in GPT-4o sometimes)
  ];
  for (const pat of patterns) {
    const m = text.match(pat);
    if (m) return m[1].trim();
  }
  return null;
}

// METRIC 1: Pass@k  — Code Executability - Unbiased estimator: pass@k = 1 - C(n-c,k)/C(n,k)
function combination(n, r) {
  if (r > n) return 0;
  if (r === 0 || r === n) return 1;
  // Use BigInt to avoid float overflow on large combos
  let res = 1n;
  for (let i = 0; i < r; i++) res = res * BigInt(n - i) / BigInt(i + 1);
  return Number(res);
}

function passAtK(n, c, k) {
  if (c === 0) return 0;
  if (n - c < k) return 1.0;
  return 1.0 - combination(n - c, k) / combination(n, k);
}

async function executeRCode(rCode, datasetDir) {
  const tmpDir = fs.mkdtempSync(path.join(CONFIG.results_dir, "rexec_"));
  const scriptPath = path.join(tmpDir, "bench.R");

  // Strip any model-generated subdirectory references since setwd already points to hg19
  rCode = rCode
    .replace(/data_dir\s*<-\s*["']path\/to\/your\/data["']/gi, 'data_dir <- "."')
    .replace(/["']path\/to\/your\/data["']/gi, '"."')
    .replace(/["']\/path\/to\/your\/data["']/gi, '"."')
    .replace(/["']your_directory["']/gi, '"."')
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']path\/to\/your\/data["']\s*\)/gi, 'Read10X(data.dir = ".")')
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']\.\/hg19\/?["']\s*\)/gi, 'Read10X(data.dir = ".")')
    .replace(/Read10X\s*\(\s*data\.dir\s*=\s*["']hg19\/?["']\s*\)/gi, 'Read10X(data.dir = ".")');

  // Mirror mainChat.js runRCode(): auto-install any library() calls before execution - allows the generated code to succeed even without pre-installed packages
  const libMatches = [...rCode.matchAll(/library\(([^)]+)\)/g)].map(m => m[1]);
  const installPrefix = libMatches.length > 0
    ? [
        `if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)`,
        ...libMatches.map(pkg =>
          `if (!requireNamespace("${pkg}", quietly=TRUE)) tryCatch(BiocManager::install("${pkg}", ask=FALSE, quiet=TRUE), error=function(e) install.packages("${pkg}", repos="https://cloud.r-project.org", quiet=TRUE))`
        )
      ].join("\n")
    : "";

  // Suppress noisy library() calls
  rCode = rCode.replace(/library\(([^)]+)\)/g, (_, pkg) => `suppressMessages(suppressWarnings(library(${pkg})))`);

  // Wrap in setwd + install prefix + sentinel output
  const fullScript = [
    `setwd("${datasetDir.replace(/\\/g, "/")}")`,
    installPrefix,
    // Suppress only messages/warnings, not errors
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
  // Keep temp dir for inspection on failure
  log(`  [R exec] "${CONFIG.r_exec}" "${scriptPath}"`);
  return new Promise((resolve) => {
    const start = Date.now();
    execFile(
      CONFIG.r_exec,
      ["--vanilla", scriptPath],
      { timeout: CONFIG.timeout_ms },
      (error, stdout, stderr) => {
        const duration_ms = Date.now() - start;
        const success = !error && stdout.includes("BENCH_EXEC_OK");
        if (!success) {
          log(`  [R exit] error=${error?.code} signal=${error?.signal} duration=${duration_ms}ms`);
          log(`  [R stderr] ${(stderr || "(empty)").slice(0, 300)}`);
          log(`  [R stdout] ${(stdout || "(empty)").slice(0, 300)}`);
        } else {
          try { fs.rmSync(tmpDir, { recursive: true }); } catch {}
        }
        resolve({ success, stderr: (stderr || "").slice(0, 400), duration_ms });
      }
    );
  });
}

async function runPassAtK(prompts, modelFn, modelName) {
  log(`\n[Pass@k] model=${modelName}  k=${CONFIG.k}  prompts=${prompts.length}`);
  const rows = [];

  for (const p of prompts) {
    const attempts = [];
    for (let i = 0; i < CONFIG.k; i++) {
      const response = await modelFn(p.prompt);
      const rCode = extractRCode(response || "");
      if (!rCode) { attempts.push({ success: false, reason: "no_code_extracted" }); continue; }
      const exec_result = await executeRCode(rCode, CONFIG.dataset_dir);
      attempts.push({ success: exec_result.success, ...exec_result });
    }
    const c = attempts.filter(a => a.success).length;
    const pass_at_k = passAtK(CONFIG.k, c, CONFIG.k);
    rows.push({ prompt_id: p.id, dataset: p.dataset, n: CONFIG.k, c, pass_at_k, attempts });
    log(`  [${p.id}] pass@${CONFIG.k}=${pass_at_k.toFixed(3)}  (${c}/${CONFIG.k} passed)`);
  }

  const avg = rows.reduce((s, r) => s + r.pass_at_k, 0) / rows.length;
  log(`  avg Pass@${CONFIG.k} = ${avg.toFixed(3)}`);
  return { rows, average: avg };
}


// METRIC 2: Workflow Success Rate - Checks for expected Seurat pipeline steps in generated code

const WORKFLOW_STEPS = {
  load:      [/Read10X|ReadMtx|read\.table|read\.csv|readRDS|Load10X_Spatial/],
  object:    [/CreateSeuratObject|as\.Seurat/],
  qc:        [/VlnPlot|PercentageFeatureSet|subset\s*\(/],
  normalize: [/NormalizeData|SCTransform/],
  features:  [/FindVariableFeatures/],
  scale:     [/ScaleData/],
  pca:       [/RunPCA/],
  neighbors: [/FindNeighbors/],
  cluster:   [/FindClusters/],
  embed:     [/RunUMAP|RunTSNE/],
};

function scoreWorkflow(rCode) {
  const steps = Object.keys(WORKFLOW_STEPS);
  if (!rCode) return { steps_found: 0, steps_total: steps.length, success_rate: 0, found: [], missing: steps };
  const found = [], missing = [];
  for (const step of steps) {
    WORKFLOW_STEPS[step].some(p => p.test(rCode)) ? found.push(step) : missing.push(step);
  }
  return { steps_found: found.length, steps_total: steps.length, success_rate: found.length / steps.length, found, missing };
}

// METRIC 3: Hallucination / Grounding Rate - Flags function calls not in the Seurat/R whitelist
const WHITELIST = new Set([
  // Seurat core
  "Read10X","Read10X_h5","ReadMtx","Load10X_Spatial",
  "CreateSeuratObject","as.Seurat",
  "PercentageFeatureSet","VlnPlot","FeatureScatter",
  "subset","NormalizeData","FindVariableFeatures","ScaleData",
  "RunPCA","ElbowPlot","DimPlot","DimHeatmap",
  "FindNeighbors","FindClusters","RunUMAP","RunTSNE",
  "FeaturePlot","RidgePlot","DotPlot","DoHeatmap",
  "FindMarkers","FindAllMarkers","VariableFeaturePlot",
  "SCTransform","PrepSCTFindMarkers",
  "SplitObject","merge","IntegrateData","FindIntegrationAnchors",
  "RunHarmony","AddMetaData","Idents","WhichCells",
  "GetAssayData","SetAssayData","Embeddings","Reductions",
  "SpatialDimPlot","SpatialFeaturePlot","FindSpatiallyVariableFeatures",
  // Base R
  "readRDS","saveRDS","read.table","read.csv","write.csv",
  "head","tail","dim","nrow","ncol","str","summary","table",
  "library","require","install.packages","BiocManager",
  "setwd","getwd","list.files","file.exists",
  "print","cat","paste","paste0","sprintf","format",
  "c","list","data.frame","matrix","as.data.frame","as.matrix","as.vector",
  "mean","sd","median","var","cor","t.test","wilcox.test",
  "grep","grepl","gsub","sub","toupper","tolower","trimws",
  "is.null","is.na","na.omit","complete.cases","which",
  "for","if","else","function","return","next","break",
  "tryCatch","try","stop","warning","message",
  "suppressMessages","suppressWarnings","suppressPackageStartupMessages",
  // ggplot2
  "ggplot","aes","geom_point","geom_bar","geom_histogram","geom_line",
  "theme","theme_bw","theme_minimal","labs","ggtitle","xlab","ylab",
  "scale_color_manual","scale_fill_manual","facet_wrap","coord_flip",
  "ggsave","pdf","png","dev.off",
]);

function scoreHallucination(rCode) {
  if (!rCode) return { grounding_rate: 0, total_calls: 0, hallucinated: [], valid: 0 };
  const pat = /\b([A-Za-z][A-Za-z0-9_.]*)\s*\(/g;
  const all = [];
  let m;
  while ((m = pat.exec(rCode)) !== null) all.push(m[1]);
  if (all.length === 0) return { grounding_rate: 1, total_calls: 0, hallucinated: [], valid: 0 };
  const hallucinated = [...new Set(all.filter(fn => !WHITELIST.has(fn)))];
  const valid = all.filter(fn => WHITELIST.has(fn)).length;
  return { grounding_rate: valid / all.length, total_calls: all.length, valid, hallucinated: hallucinated.slice(0, 20) };
}

// Run workflow + hallucination scoring for one model on all prompts
async function runCodeQualityMetrics(prompts, modelFn, modelName) {
  log(`\n[CodeQuality] model=${modelName}`);
  const workflow_rows = [], hallucination_rows = [];
  for (const p of prompts) {
    const response = await modelFn(p.prompt);
    const rCode = extractRCode(response || "");
    workflow_rows.push({ prompt_id: p.id, ...scoreWorkflow(rCode) });
    hallucination_rows.push({ prompt_id: p.id, ...scoreHallucination(rCode) });
  }
  const avg_workflow   = workflow_rows.reduce((s, r) => s + r.success_rate, 0) / workflow_rows.length;
  const avg_grounding  = hallucination_rows.reduce((s, r) => s + r.grounding_rate, 0) / hallucination_rows.length;
  log(`  avg workflow success rate: ${avg_workflow.toFixed(3)}`);
  log(`  avg grounding rate:        ${avg_grounding.toFixed(3)}`);
  return {
    workflow:     { rows: workflow_rows,     average: avg_workflow },
    hallucination:{ rows: hallucination_rows, average: avg_grounding },
  };
}

// METRIC 4: Routing Accuracy — Precision / Recall / F1
// Uses the same regex logic as orchestrator's detectDomainAndFormat()
function classifyLocally(text) {
  const t = text.toLowerCase();
  if (/\b(pubmed|literature|medical|biomedical|scientific|research|study|article|journal|doi|pmid|clinical)\b/i.test(t)) return "pubmed";
  if (/\b(r code|r script|r function|rscript|rmarkdown)\b|\b(in|using|with)\s+r\b/i.test(t)) return "code";
  if (/\b(geo|gse|dataset|database|expression|transcriptomic|sequencing|rnaseq|data|repository|accession)\b/i.test(t)) return "geo";
  if (/\b(plot|figure|chart|image|heatmap|visualization|diagram)\b/i.test(t) &&
      !/\b(r code|r script|write|generate|create)\b/i.test(t)) return "vision";
  return "general";
}

function computeF1(tp, fp, fn) {
  const p  = (tp + fp > 0) ? tp / (tp + fp) : 0;
  const r  = (tp + fn > 0) ? tp / (tp + fn) : 0;
  const f1 = (p + r > 0)   ? 2 * p * r / (p + r) : 0;
  return { precision: p, recall: r, f1 };
}

async function runRoutingAccuracy(queries) {
  log(`\n[Routing] evaluating ${queries.length} labeled queries`);
  const classes = ["pubmed", "code", "geo", "vision", "general"];
  const conf = Object.fromEntries(classes.map(c => [c, { tp: 0, fp: 0, fn: 0 }]));
  const rows = [];

  for (const q of queries) {
    const predicted = classifyLocally(q.query);
    const correct   = predicted === q.label;
    if (correct) {
      conf[q.label].tp++;
    } else {
      conf[predicted] && conf[predicted].fp++;
      conf[q.label]   && conf[q.label].fn++;
    }
    rows.push({ id: q.id, label: q.label, predicted, correct });
  }

  const perClass = Object.fromEntries(
    classes.map(c => [c, computeF1(conf[c].tp, conf[c].fp, conf[c].fn)])
  );
  const macro = {
    precision: classes.reduce((s, c) => s + perClass[c].precision, 0) / classes.length,
    recall:    classes.reduce((s, c) => s + perClass[c].recall, 0)    / classes.length,
    f1:        classes.reduce((s, c) => s + perClass[c].f1, 0)        / classes.length,
  };
  const accuracy = rows.filter(r => r.correct).length / rows.length;

  log(`  Accuracy: ${(accuracy * 100).toFixed(1)}%   Macro-F1: ${macro.f1.toFixed(3)}`);
  classes.forEach(c =>
    log(`  ${c.padEnd(8)} P=${perClass[c].precision.toFixed(2)} R=${perClass[c].recall.toFixed(2)} F1=${perClass[c].f1.toFixed(2)}`)
  );

  return { rows, perClass, macro, accuracy };
}

// METRIC 5: Genomics QA — Normalized Gene-Turing Score

function normStringMatch(pred, ref) {
  // Coerce to string safely — response may be an object in some edge cases
  const toString = v => {
    if (!v) return "";
    if (typeof v === "string") return v;
    if (typeof v === "object") return JSON.stringify(v);
    return String(v);
  };
  const norm = s => toString(s).toLowerCase().replace(/[^a-z0-9\s]/g, "").trim();
  const p = norm(pred), r = norm(ref);
  if (p === r) return 1.0;
  const pSet = new Set(p.split(/\s+/));
  const rWords = r.split(/\s+/);
  return rWords.filter(w => pSet.has(w)).length / rWords.length;
}

async function runGenomicsQA(pairs, modelFn, modelName) {
  log(`\n[GenomicsQA] model=${modelName}  pairs=${pairs.length}`);
  const rows = [];
  for (const qa of pairs) {
    const rawResponse = await modelFn(qa.question);
    // Coerce response to string
    const response = !rawResponse ? "" :
      typeof rawResponse === "string" ? rawResponse :
      typeof rawResponse === "object" && rawResponse.response ? rawResponse.response :
      JSON.stringify(rawResponse);
    const score = normStringMatch(response, qa.answer);
    rows.push({ id: qa.id, question: qa.question, expected: qa.answer, score, response: response.slice(0, 200) });
    log(`  QA ${qa.id}: score=${score.toFixed(2)}`);
  }
  const ngts = rows.reduce((s, r) => s + r.score, 0) / rows.length;
  log(`  Normalized Gene-Turing Score: ${ngts.toFixed(3)}`);
  return { rows, normalized_gene_turing_score: ngts };
}

// SUMMARY + CSV EXPORT

function printSummary() {
  const bar = "─".repeat(68);
  log(`\n${bar}\nBENCHMARK SUMMARY\n${bar}`);

  if (Object.keys(results.internal).length) {
    log("\nINTERNAL METRICS");
    log(`${"Model".padEnd(10)} ${"Pass@k".padEnd(10)} ${"Workflow".padEnd(11)} ${"Grounding".padEnd(11)}`);
    log("─".repeat(44));
    for (const [m, d] of Object.entries(results.internal)) {
      log([
        m.padEnd(10),
        `${(d.pass_at_k.average * 100).toFixed(1)}%`.padEnd(10),
        `${(d.workflow.average * 100).toFixed(1)}%`.padEnd(11),
        `${(d.hallucination.average * 100).toFixed(1)}%`.padEnd(11),
      ].join(" "));
    }
  }

  if (results.routing) {
    log(`\nROUTING ACCURACY\n  Accuracy: ${(results.routing.accuracy * 100).toFixed(1)}%  Macro-F1: ${results.routing.macro.f1.toFixed(3)}`);
  }

  if (Object.keys(results.qa).length) {
    log("\nGENOMICS QA (Gene-Turing Score)");
    for (const [m, d] of Object.entries(results.qa)) {
      log(`  ${m.padEnd(10)}: ${d.normalized_gene_turing_score.toFixed(3)}`);
    }
  }
  log(bar);
}

function writeCsv(csvPath) {
  const rows = ["model,metric,value"];
  for (const [m, d] of Object.entries(results.internal)) {
    rows.push(`${m},pass_at_k,${d.pass_at_k.average.toFixed(4)}`);
    rows.push(`${m},workflow_success,${d.workflow.average.toFixed(4)}`);
    rows.push(`${m},grounding_rate,${d.hallucination.average.toFixed(4)}`);
  }
  if (results.routing) {
    rows.push(`rolv,routing_accuracy,${results.routing.accuracy.toFixed(4)}`);
    rows.push(`rolv,routing_macro_f1,${results.routing.macro.f1.toFixed(4)}`);
  }
  for (const [m, d] of Object.entries(results.qa)) {
    rows.push(`${m},gene_turing_score,${d.normalized_gene_turing_score.toFixed(4)}`);
  }
  fs.writeFileSync(csvPath, rows.join("\n"));
  log(`CSV → ${csvPath}`);
}

// ENTRY POINT

async function main() {
  log(`${"═".repeat(60)}\nROLV Benchmark  run_id=${CONFIG.run_id}\nk=${CONFIG.k}  suite=${ARGS.suite || "all"}\n${"═".repeat(60)}`);

  const loadJsonl = (fpath) => {
    if (!fs.existsSync(fpath)) { log(`[WARN] missing: ${fpath}`); return []; }
    return fs.readFileSync(fpath, "utf-8").split("\n").filter(Boolean).map(l => JSON.parse(l));
  };

  const codePrompts = loadJsonl(CONFIG.prompts_file);
  const routingQs   = loadJsonl(CONFIG.routing_file);
  const qaPrompts   = loadJsonl(CONFIG.qa_file);
  const suite       = ARGS.suite || "all";

  // Internal suite: Pass@k + Workflow + Hallucination
  if ((suite === "all" || suite === "internal") && codePrompts.length > 0) {
    for (const [modelName, modelFn] of Object.entries(MODELS)) {
      log(`\n══ Internal Suite: ${modelName} ══`);
      const passk   = await runPassAtK(codePrompts, modelFn, modelName);
      const quality = await runCodeQualityMetrics(codePrompts, modelFn, modelName);
      results.internal[modelName] = { pass_at_k: passk, ...quality };
    }
  }

  // Routing accuracy (ROLV only — it's your classifier)
  if ((suite === "all" || suite === "routing") && routingQs.length > 0) {
    results.routing = await runRoutingAccuracy(routingQs);
  }

  // Genomics QA: all three models
  if ((suite === "all" || suite === "qa") && qaPrompts.length > 0) {
    for (const [modelName, modelFn] of Object.entries(MODELS)) {
      results.qa[modelName] = await runGenomicsQA(qaPrompts, modelFn, modelName);
    }
  }

  printSummary();
  const jsonPath = saveResults();
  writeCsv(jsonPath.replace(".json", "_summary.csv"));
  log("\n✅ Done.");
}

main().catch(err => { log(`[FATAL] ${err.stack}`); process.exit(1); });

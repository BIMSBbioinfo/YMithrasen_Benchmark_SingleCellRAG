/**
 * coding_agent_runner.js
 *
 * Automated benchmark harness for coding agents on scRNA-seq tasks.
 * Inspired by HumanEval (Pass@k), SWE-bench (repo-grounded eval),
 * Terminal-Bench (end-to-end workflow), and scBench (deterministic grading).
 *
 * Usage:
 *   node coding_agent_runner.js --agent claude
 *   node coding_agent_runner.js --agent codex
 *   node coding_agent_runner.js --agent gemini
 *   node coding_agent_runner.js --agent claude --prompts p01,p02,p03
 *   node coding_agent_runner.js --agent all
 */

import "dotenv/config";
import { exec, execSync } from "child_process";
import fs   from "fs";
import path from "path";
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
  r_exec       : process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",

  // Benchmark settings
  k            : 1,          // attempts per prompt (set to 3 for Pass@k)
  timeout_ms   : 600000,     // 10 minutes per agent call
  r_timeout_ms : 600000,     // 10 minutes per R execution
  no_recipe    : process.argv.includes("--no-recipe") || process.env.BENCH_NO_RECIPE === "1",
  context_mode : process.env.BENCH_CONTEXT_MODE || "summary", // "summary" avoids context bloat; use "full" for full scripts

  // Agent CLI commands (sandbox: recipes only + dataset read-only)
  agents: {
    claude: {
      cmd: (prompt, work_dir) =>
        `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | claude --print"`,
      extract: extractRCode,
    },
    codex: {
      cmd: (prompt, work_dir) =>
        `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | codex exec -s read-only --skip-git-repo-check --ephemeral --color never -"`,
      extract: extractRCode,
    },
    gemini: {
      // --skip-trust + default mode = read-only (write_file blocked as we confirmed)
      cmd: (prompt, work_dir) =>
        `gemini -p "${escape(prompt)}" --skip-trust --output-format text`,
      extract: extractRCode,
    },
  },
};

// ── Helpers ────────────────────────────────────────────────────────────────────

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);

function log(msg) {
  const ts = new Date().toISOString();
  process.stdout.write(`[${ts}] ${msg}\n`);
}

function escape(str) {
  return str.replace(/"/g, '\\"').replace(/\n/g, " ");
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

// ── Script preselection ────────────────────────────────────────────────────────
// Maps prompt keywords to the most relevant recipe script(s).
// Only matched scripts are injected, avoiding full recipe.md context bloat.

const SCRIPT_MAP = [
  { scripts: ["02_load_data.R"],
    keywords: ["load", "read10x", "barcodes", "matrix.mtx", "seurat object",
               "convert", "three files", "mtx", "genes.tsv", "barcodes.tsv"] },
  { scripts: ["03_qc_metrics.R"],
    keywords: ["qc", "quality control", "mitochondrial", "percent.mt",
               "violin", "nfeature", "ncount", "qc metrics", "ribosomal"] },
  { scripts: ["04_filter_cells.R"],
    keywords: ["filter", "low-quality", "low quality", "min_feat", "max_feat",
               "remove cells", "keep cells", "200 and 2500"] },
  { scripts: ["05_doublets.R"],
    keywords: ["doublet", "doubletfinder", "singlet", "multiplet"] },
  { scripts: ["06_cell_cycle.R"],
    keywords: ["cell cycle", "s phase", "g2m", "cellcyclescoring", "phase"] },
  { scripts: ["07_normalise.R"],
    keywords: ["normaliz", "normalise", "lognormalize", "sctransform", "scale.factor"] },
  { scripts: ["08_variable_genes.R"],
    keywords: ["variable", "hvg", "findvariablefeatures", "highly variable", "vst"] },
  { scripts: ["09_scale_data.R"],
    keywords: ["scale", "scaledata", "regress", "vars.to.regress"] },
  { scripts: ["10_pca.R"],
    keywords: ["pca", "elbow", "principal component", "runpca", "linear"] },
  { scripts: ["11_umap.R"],
    keywords: ["umap", "embedding", "tsne", "runumap", "visuali"] },
  { scripts: ["12_neighbours.R"],
    keywords: ["neighbour", "neighbor", "knn", "snn", "findneighbors"] },
  { scripts: ["13_clustering.R"],
    keywords: ["cluster", "resolution", "louvain", "leiden", "findclusters"] },
  { scripts: ["14_marker_genes.R"],
    keywords: ["marker", "findallmarkers", "findmarkers", "heatmap", "doheatmap",
               "differential expression", "wilcoxon", "mast", "featureplot"] },
  { scripts: ["15_cell_annotation.R"],
    keywords: ["annotate", "annotation", "cell type", "singler", "renameidents",
               "cd4", "cd8", "monocyte", "b cell", "t cell", "nk cell"] },
  { scripts: ["17_integration.R"],
    keywords: ["integrat", "harmony", "cca", "rpca", "batch correct"] },
  { scripts: ["18_pseudobulk_de.R"],
    keywords: ["pseudobulk", "deseq2", "aggregate"] },
  { scripts: ["19_trajectory.R"],
    keywords: ["trajectory", "pseudotime", "monocle"] },
];

const FULL_PIPELINE_KEYWORDS = [
  "full pipeline", "complete pipeline", "end-to-end", "full standard",
  "complete.*seurat", "entire pipeline", "from.*loading.*through",
  "complete end-to-end", "full analysis pipeline"
];

const FULL_PIPELINE_SCRIPTS = [
  "02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "05_doublets.R",
  "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R",
  "11_umap.R", "12_neighbours.R", "13_clustering.R", "14_marker_genes.R"
];

// Script dependencies (for context expansion if needed) — ensures all required scripts are included
const SCRIPT_DEPS = {
  "03_qc_metrics.R":      ["02_load_data.R"],
  "04_filter_cells.R":    ["02_load_data.R", "03_qc_metrics.R"],
  "05_doublets.R":        ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R", "11_umap.R", "12_neighbours.R", "13_clustering.R"],
  "06_cell_cycle.R":      ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R"],
  "07_normalise.R":       ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R"],
  "08_variable_genes.R":  ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R"],
  "09_scale_data.R":      ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R"],
  "10_pca.R":             ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R"],
  "11_umap.R":            ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R"],
  "12_neighbours.R":      ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R"],
  "13_clustering.R":      ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R", "12_neighbours.R", "11_umap.R"],
  "14_marker_genes.R":    ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R", "12_neighbours.R", "13_clustering.R", "11_umap.R"],
  "15_cell_annotation.R": ["02_load_data.R", "03_qc_metrics.R", "04_filter_cells.R", "07_normalise.R", "08_variable_genes.R", "09_scale_data.R", "10_pca.R", "12_neighbours.R", "13_clustering.R", "11_umap.R", "14_marker_genes.R"],
  "17_integration.R":     ["16_load_multisample.R"],
  "18_pseudobulk_de.R":   ["16_load_multisample.R", "17_integration.R", "13_clustering.R", "14_marker_genes.R", "15_cell_annotation.R"],
  "19_trajectory.R":      ["16_load_multisample.R", "17_integration.R", "13_clustering.R", "15_cell_annotation.R"],
};

// Recursively expand script dependencies to ensure all required context is included
function expandDependencies(selected) {
  const expanded = [];
  const seen = new Set();

  function add(script) {
    for (const dep of SCRIPT_DEPS[script] || []) add(dep);
    if (!seen.has(script)) {
      seen.add(script);
      expanded.push(script);
    }
  }

  selected.forEach(add);
  return expanded;
}


function loadScriptSummary(script_names) {
  return script_names
    .map(name => RECIPE_CARDS[name] || `${name} — follow the corresponding recipe script if present.`)
    .join("\n\n");
}
function selectScripts(promptText) {
  const resolved = resolveRecipeNodes(promptText);
  return resolved.scripts.length > 0 ? resolved.scripts : ["02_load_data.R"];
}

function loadScriptContent(script_names) {
  const code_dir = path.join(CONFIG.recipes_dir, "code");
  return script_names
    .map(name => {
      const fpath = path.join(code_dir, name);
      if (!fs.existsSync(fpath)) return "";
      return `# === ${name} ===\n${fs.readFileSync(fpath, "utf8")}`;
    })
    .filter(Boolean)
    .join("\n\n");
}

function extractSection(content, startHeading, endHeadingRegex) {
  const start = content.indexOf(startHeading);
  if (start === -1) return "";

  const rest = content.slice(start);
  const afterHeading = rest.slice(startHeading.length);
  const endMatch = afterHeading.match(endHeadingRegex);

  if (!endMatch || endMatch.index === undefined) {
    return rest.trim();
  }

  return rest.slice(0, startHeading.length + endMatch.index).trim();
}

function extractNotes(recipe_content) {
  const start = recipe_content.indexOf("## Notes");
  const end = recipe_content.indexOf("## Limitations");

  if (start !== -1 && end !== -1) {
    return recipe_content.slice(start, end).trim();
  }

  return "";
}

function scriptNumber(scriptName) {
  const match = scriptName.match(/^(\d+)_/);
  return match ? parseInt(match[1], 10) : null;
}

function extractRelevantSteps(recipe_content, selected_scripts) {
  const stepsSection = extractSection(
    recipe_content,
    "## Steps",
    /\n## Workflows|\n## Mid-pipeline|\n## Code references|\n## Environment|\n## Notes|\n## Limitations/
  );

  if (!stepsSection) return "";

  const selectedNums = selected_scripts
    .map(scriptNumber)
    .filter(n => Number.isInteger(n));

  if (selectedNums.length === 0) return "";

  const lines = stepsSection.split(/\r?\n/);
  const kept = [];
  let keepCurrent = false;

  for (const line of lines) {
    const stepMatch = line.match(/^\s*(\d+)\.\s+/);

    if (stepMatch) {
      const n = parseInt(stepMatch[1], 10);
      keepCurrent = selectedNums.includes(n);
    }

    if (keepCurrent) {
      kept.push(line);
    }
  }

  return kept.join("\n").trim();
}

function loadWorkflowGraph(kind = "single") {
  const filename = kind === "multi"
    ? "multi-sample-workflow.json"
    : "single-sample-workflow.json";

  const candidates = [
    path.join(CONFIG.recipes_dir, filename),
    path.join(CONFIG.workflows_dir, filename),
    path.join(CONFIG.recipes_dir, "workflows", filename),
  ];

  for (const workflowFile of candidates) {
    if (fs.existsSync(workflowFile)) {
      return JSON.parse(fs.readFileSync(workflowFile, "utf8"));
    }
  }

  log(`  [recipe] workflow file not found. Tried:`);
  for (const f of candidates) {
    log(`  [recipe]   ${f}`);
  }

  return null;
}

function detectWorkflowKind(promptText) {
  const text = promptText.toLowerCase();

  if (
    /multi[- ]?sample|multiple samples|several samples|batch|condition|donor|replicate|integrat|pseudobulk|deseq2|trajectory|pseudotime|cellchat|cell[- ]cell|pathway|enrichment/.test(text)
  ) {
    return "multi";
  }

  return "single";
}

function detectInputType(promptText, workflowKind) {
  const text = promptText.toLowerCase();

  if (/clustered object|already clustered|seurat_object_clustered/.test(text)) {
    return "seurat_object_clustered";
  }

  if (/pca object|already pca|seurat_object_pca/.test(text)) {
    return "seurat_object_pca";
  }

  if (/normalised object|normalized object|already normalised|already normalized|seurat_object_normalised|seurat_object_normalized/.test(text)) {
    return "seurat_object_normalised";
  }

  if (/qc[- ]?passed|filtered object|seurat_object_filtered/.test(text)) {
    return "seurat_object_filtered";
  }

  if (/raw seurat object|seurat_object_raw/.test(text)) {
    return "seurat_object_raw";
  }

  if (/sample list|sample_list/.test(text)) {
    return "sample_list";
  }

  if (/integrated object|seurat_object_integrated/.test(text)) {
    return "seurat_object_integrated";
  }

  if (/annotated object|seurat_object_annotated/.test(text)) {
    return "seurat_object_annotated";
  }

  if (workflowKind === "multi") {
    return "multi_sample_files";
  }

  return "raw_count_data";
}

function detectTargetNodeIds(promptText, workflowKind) {
  const text = promptText.toLowerCase();
  const targets = [];

  if (/load|read10x|matrix\.mtx|barcodes|genes\.tsv|features\.tsv|seurat object/.test(text)) {
    targets.push(workflowKind === "multi" ? "load_multisample" : "load_data");
  }

  if (/qc|quality control|mitochondrial|percent\.mt|nfeature|ncount|violin|scatter/.test(text)) {
    targets.push("qc_metrics");
  }

  if (/filter|low-quality|low quality|remove cells|keep cells|min_feat|max_feat|max_mt/.test(text)) {
    targets.push("filter_cells");
  }

  if (/doublet|doubletfinder|singlet|multiplet/.test(text)) {
    targets.push("doublets");
  }

  if (/cell cycle|s phase|g2m|phase/.test(text)) {
    targets.push("cell_cycle");
  }

  if (/normaliz|normalise|lognormalize|sctransform/.test(text)) {
    targets.push("normalise");
  }

  if (/variable|hvg|highly variable|findvariablefeatures/.test(text)) {
    targets.push("variable_genes");
  }

  if (/scale|scaledata|regress/.test(text)) {
    targets.push("scale_data");
  }

  if (/pca|principal component|elbow|runpca/.test(text)) {
    targets.push("pca");
  }

  if (/umap|embedding|tsne|visuali|dimplot/.test(text)) {
    targets.push("umap");
  }

  if (/neighbo[u]?r|knn|snn|findneighbors/.test(text)) {
    targets.push("neighbours");
  }

  if (/cluster|resolution|louvain|leiden|findclusters/.test(text)) {
    targets.push("clustering");
  }

  if (/marker|findallmarkers|findmarkers|differential expression|wilcoxon|mast/.test(text)) {
    targets.push("marker_genes");
  }

  if (/annotat|cell type|singler|renameidents/.test(text)) {
    targets.push("cell_annotation");
  }

  if (/integrat|batch|harmony|cca|rpca/.test(text)) {
    targets.push("integration");
  }

  if (/pseudobulk|deseq2|aggregate/.test(text)) {
    targets.push("pseudobulk_de");
  }

  if (/trajectory|pseudotime|monocle/.test(text)) {
    targets.push("trajectory");
  }

  if (/cellchat|cell[- ]cell communication|communication network/.test(text)) {
    targets.push("cellchat");
  }

  if (/pathway|enrichr|go|kegg|msigdb|enrichment/.test(text)) {
    targets.push("pathway_enrichment");
  }

  if (targets.length === 0) {
    targets.push(workflowKind === "multi" ? "load_multisample" : "load_data");
  }

  return [...new Set(targets)];
}

function buildAdjacency(workflow) {
  const adj = new Map();

  for (const node of workflow.nodes || []) {
    adj.set(node.id, []);
  }

  for (const edge of workflow.edges || []) {
    if (!adj.has(edge.from)) adj.set(edge.from, []);
    adj.get(edge.from).push(edge.to);
  }

  return adj;
}

function nodeForInputType(workflow, inputType) {
  const exact = (workflow.nodes || []).find(n => n.input_type === inputType);
  if (exact) return exact.id;

  const producedBy = (workflow.nodes || []).find(n => n.output_type === inputType);
  if (producedBy) return producedBy.id;

  return null;
}

function shortestPath(workflow, startNodeId, targetNodeId) {
  if (!startNodeId || !targetNodeId) return [];
  if (startNodeId === targetNodeId) return [targetNodeId];

  const adj = buildAdjacency(workflow);
  const queue = [[startNodeId]];
  const seen = new Set([startNodeId]);

  while (queue.length > 0) {
    const pathSoFar = queue.shift();
    const last = pathSoFar[pathSoFar.length - 1];

    for (const next of adj.get(last) || []) {
      if (seen.has(next)) continue;

      const nextPath = [...pathSoFar, next];

      if (next === targetNodeId) {
        return nextPath;
      }

      seen.add(next);
      queue.push(nextPath);
    }
  }

  return [targetNodeId];
}

function resolveRecipeNodes(promptText) {
  // Full pipeline shortcut — if the prompt matches a full pipeline pattern,
  // return all standard single-sample scripts immediately without keyword matching
  if (FULL_PIPELINE_KEYWORDS.some(kw => new RegExp(kw, "i").test(promptText))) {
    return {
      workflowKind : "single",
      inputType    : "raw_count_data",
      startNode    : "load_data",
      targetNodes  : ["marker_genes"],
      nodes        : [],
      scripts      : FULL_PIPELINE_SCRIPTS,
    };
  }

  const workflowKind = detectWorkflowKind(promptText);
  const workflow = loadWorkflowGraph(workflowKind);

  if (!workflow) {
    const fallbackTargetNodes = detectTargetNodeIds(promptText, workflowKind);

    return {
      workflowKind,
      inputType: workflowKind === "multi" ? "multi_sample_files" : "raw_count_data",
      startNode: workflowKind === "multi" ? "load_multisample" : "load_data",
      targetNodes: fallbackTargetNodes,
      nodes: [],
      scripts: workflowKind === "multi" ? ["16_load_multisample.R"] : ["02_load_data.R"],
    };
  }

  const inputType = detectInputType(promptText, workflowKind);
  const startNode = nodeForInputType(workflow, inputType);
  const targetNodes = detectTargetNodeIds(promptText, workflowKind);

  const nodeIds = [];

  for (const target of targetNodes) {
    const pathNodes = shortestPath(workflow, startNode, target);

    for (const nodeId of pathNodes) {
      if (!nodeIds.includes(nodeId)) {
        nodeIds.push(nodeId);
      }
    }
  }

  const selectedNodes = nodeIds
    .map(id => (workflow.nodes || []).find(n => n.id === id))
    .filter(Boolean);

  const scripts = selectedNodes
    .map(n => n.path)
    .filter(Boolean)
    .map(p => path.basename(p));

  return {
    workflowKind,
    inputType,
    startNode,
    targetNodes,
    nodes: selectedNodes,
    scripts: [...new Set(scripts)],
  };
}

function extractParamsTable(recipe_content) {
  const start = recipe_content.indexOf("| Parameter |");
  const end = recipe_content.indexOf("## Steps");

  if (start !== -1 && end !== -1) {
    return recipe_content.slice(start, end).trim();
  }

  return "";
}

function extractNotesForSelectedNodes(recipe_content, selectedNodes) {
  const selectedIds = new Set(selectedNodes.map(n => n.id));
  const notes = [];

  if (selectedIds.has("qc_metrics") || selectedIds.has("filter_cells")) {
    notes.push("- For human data, use `^MT-` for mitochondrial genes. Compute `percent.mt` before filtering by mitochondrial percentage.");
  }

  if (selectedIds.has("doublets")) {
    notes.push("- Doublet detection is optional-package-dependent. Use DoubletFinder only if it is installed. Do not install packages inside the benchmark script.");
    notes.push("- If DoubletFinder is unavailable, create a safe fallback output that marks all cells as `Singlet` and continues, rather than failing.");
    notes.push("- If using DoubletFinder, detect classification columns dynamically with `grep()`, not hardcoded names.");
  }

  if (selectedIds.has("cell_annotation")) {
    notes.push("- Prefer manual marker-based annotation for benchmark robustness.");
    notes.push("- Use SingleR/celldex only if explicitly requested and installed. Do not install packages inside the benchmark script.");
    notes.push("- If no reference package is available, assign simple cluster labels or marker-informed labels and save the annotated object/output required by the task.");
  }

  if (selectedIds.has("marker_genes")) {
    notes.push("- Marker genes require clusters to exist first. Use FindAllMarkers only after FindClusters has run.");
  }

  if (selectedIds.has("integration")) {
    notes.push("- Multi-sample integration should use common genes across samples and avoid optional Harmony unless explicitly requested and installed.");
  }

  if (selectedIds.has("trajectory")) {
    notes.push("- Use Monocle3 only if explicitly requested and installed. Otherwise avoid trajectory code that will fail.");
  }

  if (selectedIds.has("cellchat")) {
    notes.push("- Use CellChat only if explicitly requested and installed. Otherwise avoid CellChat code that will fail.");
  }

  return notes.join("\n");
}

function buildWorkflowNodeContext(resolved) {
  if (!resolved.nodes || resolved.nodes.length === 0) {
    return "No workflow nodes selected.";
  }

  return resolved.nodes.map((node, idx) => {
    return `${idx + 1}. ${node.id} (${node.input_type} -> ${node.output_type})
   script: ${node.path}`;
  }).join("\n");
}

// Pointer-style cards for each recipe script.
// Each card tells the agent: what the script does, which function to use, and where to find it.
// No prose summaries — the actual file is the source of truth.
const RECIPE_POINTERS = {
  "02_load_data.R": `## 02_load_data.R — Load raw count data
File: ./recipes/code/02_load_data.R
Use Read10X() for MTX folders, Read10X_h5() for H5 files, readRDS() for RDS.
Creates a Seurat object with CreateSeuratObject(). Object named: seu`,

  "03_qc_metrics.R": `## 03_qc_metrics.R — Compute QC metrics
File: ./recipes/code/03_qc_metrics.R
Use PercentageFeatureSet(seu, pattern="^MT-") for percent.mt.
Adds nFeature_RNA, nCount_RNA, percent.mt to seu@meta.data.`,

  "04_filter_cells.R": `## 04_filter_cells.R — Filter low-quality cells
File: ./recipes/code/04_filter_cells.R
Use subset() with nFeature_RNA, nCount_RNA, percent.mt thresholds.
Output object named: seu_filtered`,

  "05_doublets.R": `## 05_doublets.R — Doublet detection
File: ./recipes/code/05_doublets.R
Key function: run_doubletfinder_safe(seu, n_pcs, doublet_rate, pK, run_pk_sweep)
Usage: source("./recipes/code/05_doublets.R") then call run_doubletfinder_safe()
DO NOT call doubletFinder() directly — use run_doubletfinder_safe() only.
Input: seu_filtered with clusters and UMAP already run.
Output: seu_df with $doublet_status; subset to seu_singlets after.`,

  "06_cell_cycle.R": `## 06_cell_cycle.R — Cell cycle scoring
File: ./recipes/code/06_cell_cycle.R
Use CellCycleScoring() with s.genes and g2m.genes from Seurat.
Adds S.Score, G2M.Score, Phase to seu@meta.data.`,

  "07_normalise.R": `## 07_normalise.R — Normalise data
File: ./recipes/code/07_normalise.R
Use NormalizeData(seu_filtered, normalization.method="LogNormalize", scale.factor=10000).
Output object named: seu_filtered (updated in place).`,

  "08_variable_genes.R": `## 08_variable_genes.R — Find highly variable genes
File: ./recipes/code/08_variable_genes.R
Use FindVariableFeatures(seu_filtered, selection.method="vst", nfeatures=2000).`,

  "09_scale_data.R": `## 09_scale_data.R — Scale data
File: ./recipes/code/09_scale_data.R
Use ScaleData(seu_filtered, features=VariableFeatures(seu_filtered)).`,

  "10_pca.R": `## 10_pca.R — PCA
File: ./recipes/code/10_pca.R
Use RunPCA(seu_filtered, features=VariableFeatures(seu_filtered), npcs=50).
Use ElbowPlot() to choose n_pcs.`,

  "11_umap.R": `## 11_umap.R — UMAP embedding
File: ./recipes/code/11_umap.R
Use RunUMAP(seu_filtered, dims=1:n_pcs).`,

  "12_neighbours.R": `## 12_neighbours.R — Find neighbours
File: ./recipes/code/12_neighbours.R
Use FindNeighbors(seu_filtered, dims=1:n_pcs).`,

  "13_clustering.R": `## 13_clustering.R — Cluster cells
File: ./recipes/code/13_clustering.R
Use FindClusters(seu_filtered, resolution=0.5).
Cluster labels available as Idents(seu_filtered) or seu_filtered$seurat_clusters.`,

  "14_marker_genes.R": `## 14_marker_genes.R — Find marker genes
File: ./recipes/code/14_marker_genes.R
Use FindAllMarkers(seu_filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25).
Requires clusters to exist first (FindClusters must have run).`,

  "15_cell_annotation.R": `## 15_cell_annotation.R — Annotate cell types
File: ./recipes/code/15_cell_annotation.R
Prefer marker-based annotation using RenameIdents().
Use SingleR/celldex only if explicitly requested and installed.`,

  "16_load_multisample.R": `## 16_load_multisample.R — Load multiple samples
File: ./recipes/code/16_load_multisample.R
Read each sample with Read10X(), create Seurat objects, merge with merge().`,

  "17_integration.R": `## 17_integration.R — Integrate samples
File: ./recipes/code/17_integration.R
Use IntegrateLayers() or Harmony. Requires merged object with batch info.`,

  "18_pseudobulk_de.R": `## 18_pseudobulk_de.R — Pseudobulk differential expression
File: ./recipes/code/18_pseudobulk_de.R
Aggregate counts per sample/cluster, then use DESeq2.`,

  "19_trajectory.R": `## 19_trajectory.R — Trajectory analysis
File: ./recipes/code/19_trajectory.R
Use Monocle3 only if explicitly requested and installed.`,
};

function buildPointerCards(selected_scripts, resolved) {
  const paramsTable = (() => {
    const recipe_content = fs.existsSync(CONFIG.recipe_file)
      ? fs.readFileSync(CONFIG.recipe_file, "utf8") : "";
    return extractParamsTable(recipe_content);
  })();

  const cards = selected_scripts
    .map(name => RECIPE_POINTERS[name] || `## ${name}\nFile: ./recipes/code/${name}`)
    .join("\n\n");

  const selectedNotes = extractNotesForSelectedNodes("", resolved.nodes);

  return `RECIPE / SKILL ROUTING CONTEXT

Workflow selected: ${resolved.workflowKind}
Detected input type: ${resolved.inputType || "unknown"}
Target node(s): ${(resolved.targetNodes || []).join(", ") || "unknown"}

The selected workflow path below is dependency guidance, not mandatory code to copy.
Use it to understand what must exist before the requested target step.
Generate the smallest complete executable R script for the benchmark task.

PARAMETER GUIDANCE
The recipe parameter table provides suggested ranges only.
Use task-specified values first. Otherwise choose reasonable values and define them in code.

${paramsTable}

SELECTED WORKFLOW PATH
${selected_scripts.map((s, i) => `${i + 1}. ${s}`).join("\n")}

RECIPE CARDS — read the listed file before writing code for that step:
${cards}
${selectedNotes ? `\nTASK-SPECIFIC NOTES\n${selectedNotes}` : ""}`;
}

function buildPromptWithContext(prompt_obj) {
  if (CONFIG.no_recipe) {
    return `You are an expert bioinformatician working with single-cell RNA-seq data in R using Seurat.

DATASET PATH: ${CONFIG.dataset_dir}

TASK: ${prompt_obj.prompt}

INSTRUCTIONS:
- Use appropriate Seurat functions for scRNA-seq analysis
- Set seed = 42 for reproducibility
- Save any output files to the current working directory
- Write complete, executable R code
- Output ONLY a single \`\`\`r code block containing all the R code
- Do NOT ask for permission to write files
- Do NOT use tools or execute anything
- Just output the R code block and nothing else`;
  }

  const resolved = resolveRecipeNodes(prompt_obj.prompt);
  const selected_scripts = resolved.scripts;

  const recipe_context = CONFIG.context_mode === "full"
    ? loadScriptContent(selected_scripts)
    : buildPointerCards(selected_scripts, resolved);

  log(`  [recipe] workflow=${resolved.workflowKind || "unknown"}`);
  log(`  [recipe] input=${resolved.inputType || "unknown"}`);
  log(`  [recipe] targets=${(resolved.targetNodes || []).join(", ") || "unknown"}`);
  log(`  [recipe] selected: ${(selected_scripts || []).join(", ") || "none"}`);
  log(`  [recipe] mode=${CONFIG.context_mode}`);

  return `You are an expert bioinformatician working with single-cell RNA-seq data in R using Seurat.

A recipe/skill package is available in the workspace at:

./recipes/

Before writing code:
- Read ./recipes/recipe.md.
- Use the selected workflow path below to identify the relevant recipe step(s).
- For tasks matching a specific recipe step, inspect the corresponding file in ./recipes/code/ before writing code.
- Use ./recipes/single-sample-workflow.json or ./recipes/multi-sample-workflow.json only if needed.
- Treat recipe parameters as guidance, not hardcoded values.
- Use task-specified values first.
- Generate the smallest complete executable R script that satisfies the task.

SELECTED SKILL CONTEXT:
${recipe_context}

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
- Save required outputs to the current working directory.`;
}

// ── Working directory setup ────────────────────────────────────────────────────

function setupWorkDir(agent, prompt_id, attempt) {
  const work_dir = path.join(CONFIG.work_base, agent, `${prompt_id}_attempt${attempt}`);
  fs.mkdirSync(work_dir, { recursive: true });

  // Copy recipes folder (read access only — agents read from here)
  const recipes_dest = path.join(work_dir, "recipes");
  if (!fs.existsSync(recipes_dest)) {
    fs.cpSync(CONFIG.recipes_dir, recipes_dest, { recursive: true });
  }

  return work_dir;
}

function cleanWorkDir(work_dir) {
  // Remove R scripts and outputs but keep recipes
  const files = fs.readdirSync(work_dir);
  for (const f of files) {
    if (f !== "recipes") {
      fs.rmSync(path.join(work_dir, f), { recursive: true, force: true });
    }
  }
}

// ── Agent call ─────────────────────────────────────────────────────────────────

function callAgent(agent_name, prompt_text, work_dir) {
  return new Promise((resolve) => {
    const agent = CONFIG.agents[agent_name];
    const cmd = agent.cmd(prompt_text, work_dir);
    const env = {
      ...process.env,
      GEMINI_CLI_TRUST_WORKSPACE: "true",
    };

    fs.writeFileSync(path.join(work_dir, "final_prompt.txt"), prompt_text, "utf8");

    log(`  [prompt] chars=${prompt_text.length}`);
    log(`  [${agent_name}] calling agent...`);

    exec(cmd, {
      timeout : CONFIG.timeout_ms,
      cwd     : work_dir,
      env,
      input   : prompt_text,
      maxBuffer: 20 * 1024 * 1024,
    }, (error, stdout, stderr) => {
      fs.writeFileSync(path.join(work_dir, "agent_stdout.txt"), stdout || "", "utf8");
      fs.writeFileSync(path.join(work_dir, "agent_stderr.txt"), stderr || "", "utf8");

      if (error) {
        log(`  [${agent_name}] agent error: ${error.message?.slice(0, 180)}`);
      }

      resolve({
        stdout: stdout || "",
        stderr: stderr || "",
        error: error?.message || null
      });
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
          grounding_rate : 0.0,
          grounding_fns  : [],
        });
      }
    });
  });
}

// ── Single prompt evaluation ───────────────────────────────────────────────────

async function evaluatePrompt(agent_name, prompt_obj, attempt = 1) {
  const { id: prompt_id } = prompt_obj;
  log(`\n  [${agent_name}/${prompt_id}] attempt ${attempt}`);

  const work_dir    = setupWorkDir(agent_name, prompt_id, attempt);
  const prompt_text = buildPromptWithContext(prompt_obj);
  const t_start     = Date.now();

  // 1. Call agent
  const agent_result = await callAgent(agent_name, prompt_text, work_dir);
  const t_agent      = Date.now() - t_start;
  log(`  [${agent_name}/${prompt_id}] agent responded in ${(t_agent/1000).toFixed(1)}s`);

  // Handle agent CLI errors (no output + error message)
  if (agent_result.error && !agent_result.stdout) {
    return {
      agent: agent_name,
      prompt_id,
      category: prompt_obj.category,
      difficulty: prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : Date.now() - t_start,
      r_executed       : false,
      r_success        : false,
      grader_pass      : false,
      grader_score     : 0.0,
      grounding_rate   : 0.0,
      failure_type     : "agent_cli_error",
      reason           : agent_result.error,
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
      category: prompt_obj.category,
      difficulty: prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : Date.now() - t_start,
      r_executed       : false,
      r_success        : false,
      grader_pass      : false,
      grader_score     : 0.0,
      grounding_rate   : 0.0,
      failure_type     : "no_code_extracted",
      reason           : "No R code in agent output. Inspect agent_stdout.txt and agent_stderr.txt in the workdir.",
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
      category       : prompt_obj.category,
      difficulty     : prompt_obj.difficulty,
      attempt,
      agent_latency_ms : t_agent,
      total_ms         : t_total,
      r_executed       : true,
      r_success        : false,
      r_error          : r_result.error,
      grader_pass      : false,
      grader_score     : 0.0,
      grounding_rate   : 0.0,
      failure_type     : r_result.error && /timed out|timeout/i.test(r_result.error) ? "timeout" : "r_runtime_error",
      reason           : r_result.error || "R execution failed. Inspect *_r_stdout.txt in the workdir.",
      work_dir,
    };
  }
  // 4. Run grader
  const grader_result = await runGrader(prompt_id, work_dir);
  log(`  [${agent_name}/${prompt_id}] grader: pass=${grader_result.pass} score=${grader_result.score?.toFixed(2)} grounding=${grader_result.grounding_rate?.toFixed(2)}`);

  return {
    agent          : agent_name,
    prompt_id,
    category       : prompt_obj.category,
    difficulty     : prompt_obj.difficulty,
    attempt,
    agent_latency_ms : t_agent,
    total_ms         : t_total,
    r_executed       : true,
    r_success        : r_result.success,
    r_error          : r_result.error,
    grader_pass      : grader_result.pass   || false,
    grader_score     : grader_result.score  || 0.0,
    grounding_rate   : grader_result.grounding_rate || 0.0,
    grounding_fns    : grader_result.grounding_fns  || [],
    failure_type     : grader_result.pass ? "pass" : "grader_fail",
    reason           : grader_result.reason || "",
    work_dir,
  };
}

// ── Pass@k calculation ─────────────────────────────────────────────────────────

function computePassAtK(attempts) {
  // Pass@k: 1 if any attempt passes, 0 otherwise
  return attempts.some(a => a.grader_pass) ? 1.0 : 0.0;
}

// ── Main benchmark runner ──────────────────────────────────────────────────────

async function runBenchmark(agent_name, prompt_ids = null) {
  log(`\n${"═".repeat(70)}`);
  log(`Coding Agent Benchmark  run_id=${run_id}`);
  log(`Agent: ${agent_name}  k=${CONFIG.k}`);
  log(`${"═".repeat(70)}`);

  const prompts = loadJSONL(CONFIG.prompts_file);
  const selected = prompt_ids
    ? prompts.filter(p => prompt_ids.includes(p.id))
    : prompts;

  log(`Loaded ${selected.length} prompts`);
  fs.mkdirSync(CONFIG.results_dir, { recursive: true });
  fs.mkdirSync(CONFIG.work_base,   { recursive: true });

  const all_results = [];

  for (const prompt_obj of selected) {
    const attempts = [];
    for (let k = 1; k <= CONFIG.k; k++) {
      const result = await evaluatePrompt(agent_name, prompt_obj, k);
      attempts.push(result);
      all_results.push(result);
      if (result.grader_pass) break; // Stop early if passed
    }

    const pass_at_k = computePassAtK(attempts);
    log(`  [${prompt_obj.id}] Pass@${CONFIG.k}=${pass_at_k.toFixed(3)}`);
  }

  // ── Summary ──────────────────────────────────────────────────────────────────
  log(`\n${"─".repeat(70)}`);
  log(`BENCHMARK SUMMARY — ${agent_name}`);
  log(`${"─".repeat(70)}`);

  // By metric
  const pass_rate      = all_results.filter(r => r.grader_pass).length / all_results.length;
  const r_success_rate = all_results.filter(r => r.r_success).length / all_results.length;
  const avg_grounding  = all_results.reduce((s, r) => s + (r.grounding_rate || 0), 0) / all_results.length;
  const avg_latency    = all_results.reduce((s, r) => s + (r.agent_latency_ms || 0), 0) / all_results.length;

  log(`Pass@k (grader):      ${(pass_rate * 100).toFixed(1)}%`);
  log(`Workflow success (R):  ${(r_success_rate * 100).toFixed(1)}%`);
  log(`Avg grounding rate:    ${(avg_grounding * 100).toFixed(1)}%`);
  log(`Avg agent latency:     ${(avg_latency / 1000).toFixed(1)}s`);

  // By category
  const categories = [...new Set(all_results.map(r => r.category).filter(Boolean))];
  if (categories.length > 0) {
    log(`\nBy category:`);
    for (const cat of categories) {
      const cat_results = all_results.filter(r => r.category === cat);
      const cat_pass    = cat_results.filter(r => r.grader_pass).length / cat_results.length;
      log(`  ${(cat || "unknown").padEnd(20)} ${(cat_pass * 100).toFixed(1)}% (${cat_results.length} prompts)`);
    }
  }

  // By difficulty
  log(`\nBy difficulty:`);
  for (const diff of ["simple", "medium", "hard"]) {
    const diff_results = all_results.filter(r => (r.difficulty === diff) ||
      (r.prompt_id && r.prompt_id.startsWith(diff.slice(0,3))));
    if (diff_results.length === 0) continue;
    const diff_pass = diff_results.filter(r => r.grader_pass).length / diff_results.length;
    log(`  ${diff.padEnd(10)} ${(diff_pass * 100).toFixed(1)}% (${diff_results.length} prompts)`);
  }

  log(`\nBy failure type:`);

  const failure_types = [
    ...new Set(
      all_results.map(r => r.failure_type || (r.grader_pass ? "pass" : "unknown"))
    )
  ];

  for (const ft of failure_types) {
    const n = all_results.filter(r =>
      (r.failure_type || (r.grader_pass ? "pass" : "unknown")) === ft
    ).length;

    log(`  ${ft.padEnd(20)} ${n}`);
  }

  // Save results
  const out_file = path.join(CONFIG.results_dir,
    `coding_agent_${agent_name}_${run_id}.json`);
  const summary = {
    run_id, agent: agent_name,
    timestamp   : new Date().toISOString(),
    n_prompts   : selected.length,
    k           : CONFIG.k,
    metrics     : { pass_rate, r_success_rate, avg_grounding, avg_latency_ms: avg_latency },
    by_category : Object.fromEntries(categories.map(cat => {
      const r = all_results.filter(x => x.category === cat);
      return [cat, r.filter(x => x.grader_pass).length / r.length];
    })),
    results: all_results,
  };
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
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
  dataset_dir  : process.env.BENCH_DATASET_DIR ||
                 "D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19",
  prompts_file : path.join(__dirname, "prompts", "coding_agent_prompts.jsonl"),
  grader_script: path.join(__dirname, "graders", "graders.R"),
  results_dir  : path.join(__dirname, "results"),
  work_base    : path.join(__dirname, "workdirs"),
  r_exec       : process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",

  // Benchmark settings
  k            : 1,          // attempts per prompt (set to 3 for Pass@k)
  parallel     : 1,          // number of prompts to run concurrently (set via --parallel)
  no_recipe    : false,      // set true to run without recipe context (baseline condition)
  timeout_ms   : 600000,     // 10 minutes per agent call
  r_timeout_ms : 600000,     // 10 minutes per R execution

  // Agent CLI commands (sandbox: recipes only + dataset read-only)
  agents: {
    claude: {
      // --allowedTools restricts Claude Code to read-only file access
      cmd: (prompt, work_dir) =>
        `claude --print "${escape(prompt)}"`,
      extract: extractRCode,
    },
    codex: {
      // codex exec defaults to read-only sandbox
      cmd: (prompt, work_dir) =>
        `codex exec "${escape(prompt)}"`,
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

function extractRCode(text) {
  // Extract R code blocks from agent output
  // Try ```r ... ``` first, then ```{r} ... ```, then raw R patterns
  const patterns = [
    /```r\n([\s\S]*?)```/gi,
    /```\{r\}\n([\s\S]*?)```/gi,
    /```R\n([\s\S]*?)```/gi,
    /```\n(library[\s\S]*?)```/gi,
  ];
  for (const pattern of patterns) {
    const matches = [...text.matchAll(pattern)];
    if (matches.length > 0) {
      return matches.map(m => m[1]).join("\n\n");
    }
  }
  // Fallback: if output looks like R code directly (starts with library or #)
  const lines = text.split("\n");
  const r_lines = lines.filter(l =>
    l.match(/^library\(|^#|^seu|^data|^set\.seed|^cat\(|^print\(|^Read10X/)
  );
  if (r_lines.length > 5) return text;
  return null;
}

function buildPromptWithContext(prompt_obj) {
  if (CONFIG.no_recipe) {
    // No-recipe baseline condition: raw prompt only, no recipe guidance
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

  // Recipe-grounded condition: inject recipe.md as context
  const recipe_content = fs.existsSync(path.join(CONFIG.recipes_dir, "recipe.md"))
    ? fs.readFileSync(path.join(CONFIG.recipes_dir, "recipe.md"), "utf8")
    : "";

  return `You are an expert bioinformatician. You have access to an scRNA-seq analysis recipe.

RECIPE CONTEXT:
${recipe_content.slice(0, 3000)}

DATASET PATH: ${CONFIG.dataset_dir}

TASK: ${prompt_obj.prompt}

INSTRUCTIONS:
- Follow the recipe structure and parameter block pattern
- Use the exact functions specified in the recipe
- Set seed = 42 for reproducibility
- Save any output files to the current working directory
- Write complete, executable R code
- Output ONLY a single \`\`\`r code block containing all the R code
- Do NOT ask for permission to write files
- Do NOT use tools or execute anything
- Do NOT explain the code before or after
- Just output the R code block and nothing else`;
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
    const agent   = CONFIG.agents[agent_name];
    const cmd     = agent.cmd(prompt_text, work_dir);
    const env     = {
      ...process.env,
      // Gemini trust
      GEMINI_CLI_TRUST_WORKSPACE: "true",
    };

    log(`  [${agent_name}] calling agent...`);

    exec(cmd, {
      timeout : CONFIG.timeout_ms,
      cwd     : work_dir,
      env,
      maxBuffer: 10 * 1024 * 1024,
    }, (error, stdout, stderr) => {
      if (error && !stdout) {
        log(`  [${agent_name}] agent error: ${error.message?.slice(0, 100)}`);
        resolve({ stdout: "", stderr: stderr || "", error: error.message });
        return;
      }
      resolve({ stdout: stdout || "", stderr: stderr || "", error: null });
    });
  });
}

// ── R execution ───────────────────────────────────────────────────────────────

function executeRCode(r_code, work_dir, prompt_id) {
  return new Promise((resolve) => {
    const r_file    = path.join(work_dir, `${prompt_id}_generated.R`);
    const log_file  = path.join(work_dir, `${prompt_id}_r_output.txt`);
    fs.writeFileSync(r_file, r_code, "utf8");

    log(`  [R] executing ${path.basename(r_file)}`);

    exec(
      `"${CONFIG.r_exec}" "${r_file}" 2>&1`,
      { timeout: CONFIG.r_timeout_ms, cwd: work_dir, maxBuffer: 10 * 1024 * 1024 },
      (error, stdout) => {
        const combined = stdout || "";
        // Save R output to file so grader can inspect it
        try { fs.writeFileSync(log_file, combined, "utf8"); } catch {}
        const success = !error || error.code === 0;
        if (!success) {
          // Log first error line for immediate visibility
          const err_line = combined.split("\n").find(l => l.includes("Error") || l.includes("error")) || "";
          if (err_line) log(`  [R error] ${err_line.slice(0, 120)}`);
        }
        resolve({
          success,
          stdout : combined,
          error  : error?.message || null,
          r_file,
        });
      }
    );
  });
}

// ── Grader ─────────────────────────────────────────────────────────────────────

function runGrader(prompt_id, work_dir) {
  return new Promise((resolve) => {
    // Normalize to forward slashes for R compatibility on Windows
    const work_dir_r = work_dir.replace(/\\/g, "/");
    const cmd = `"${CONFIG.r_exec}" "${CONFIG.grader_script}" "${prompt_id}" "${work_dir_r}"`;
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

  // 2. Extract R code
  const r_code = extractRCode(agent_result.stdout);
  if (!r_code) {
    log(`  [${agent_name}/${prompt_id}] no R code extracted`);
    return {
      agent: agent_name, prompt_id, attempt,
      agent_latency_ms : t_agent,
      r_executed       : false,
      r_success        : false,
      grader_pass      : false,
      grader_score     : 0.0,
      grounding_rate   : 0.0,
      reason           : "No R code in agent output",
    };
  }

  // 3. Execute R code
  const r_result = await executeRCode(r_code, work_dir, prompt_id);
  const t_total  = Date.now() - t_start;
  log(`  [${agent_name}/${prompt_id}] R ${r_result.success ? "✓" : "✗"} (${(t_total/1000).toFixed(1)}s)`);

  // 4. Run grader
  const grader_result = await runGrader(prompt_id, work_dir);
  log(`  [${agent_name}/${prompt_id}] grader: pass=${grader_result.pass} score=${grader_result.score?.toFixed(2)} grounding=${grader_result.grounding_rate?.toFixed(2)}`);

  return {
    agent            : agent_name,
    prompt_id,
    category         : prompt_obj.category,
    difficulty       : prompt_obj.difficulty,
    attempt,
    agent_latency_ms : t_agent,
    total_ms         : t_total,
    r_executed       : true,
    r_success        : r_result.success,
    r_error          : r_result.error,
    r_stdout_preview : (r_result.stdout || "").slice(0, 300),
    grader_pass      : grader_result.pass   || false,
    grader_score     : grader_result.score  || 0.0,
    grounding_rate   : grader_result.grounding_rate || 0.0,
    grounding_fns    : grader_result.grounding_fns  || [],
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

  // Run prompts in parallel batches of CONFIG.parallel
  const chunks = [];
  for (let i = 0; i < selected.length; i += CONFIG.parallel) {
    chunks.push(selected.slice(i, i + CONFIG.parallel));
  }

  for (const chunk of chunks) {
    // Run all prompts in this chunk concurrently
    const chunk_results = await Promise.all(chunk.map(async (prompt_obj) => {
      const attempts = [];
      for (let k = 1; k <= CONFIG.k; k++) {
        const result = await evaluatePrompt(agent_name, prompt_obj, k);
        attempts.push(result);
        if (result.grader_pass) break; // Stop early if passed
      }
      const pass_at_k = computePassAtK(attempts);
      log(`  [${prompt_obj.id}] Pass@${CONFIG.k}=${pass_at_k.toFixed(3)}`);
      return attempts;
    }));
    // Flatten attempts into all_results
    chunk_results.forEach(attempts => attempts.forEach(r => all_results.push(r)));
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

  // Save results
  const mode_suffix = CONFIG.no_recipe ? "_no-recipe" : "_with-recipe";
  const out_file = path.join(CONFIG.results_dir,
    `coding_agent_${agent_name}${mode_suffix}_${run_id}.json`);
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

// --no-recipe flag: run without recipe context (baseline condition)
CONFIG.no_recipe = raw_args.includes("--no-recipe") || process.env.BENCH_NO_RECIPE === "true";
if (CONFIG.no_recipe) log("Mode: NO-RECIPE baseline (recipe context disabled)");

// --k flag: override k value
const k_idx = raw_args.indexOf("--k");
if (k_idx >= 0) CONFIG.k = parseInt(raw_args[k_idx + 1]) || 1;

// --parallel flag: number of concurrent prompts
const par_idx = raw_args.indexOf("--parallel");
if (par_idx >= 0) CONFIG.parallel = parseInt(raw_args[par_idx + 1]) || 1;
if (CONFIG.parallel > 1) log(`Parallel mode: ${CONFIG.parallel} concurrent prompts`);

if (agent_arg === "all") {
  for (const agent of ["claude", "codex", "gemini"]) {
    await runBenchmark(agent, prompt_ids);
  }
} else {
  await runBenchmark(agent_arg, prompt_ids);
}
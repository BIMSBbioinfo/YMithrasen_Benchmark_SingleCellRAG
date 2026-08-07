/**
 * coding_agent_skills_runner.js
 *
 * Native-skills benchmark for coding agents on scRNA-seq tasks.
 *
 * Experimental condition: oracle-activated native progressive-disclosure skill.
 *   - Claude Code: .claude/skills/single-cell-r/SKILL.md  (invoked via /single-cell-r)
 *   - Codex:       .agents/skills/single-cell-r/SKILL.md  (invoked via explicit name)
 *   - Full recipe body NOT injected into the task prompt
 *   - Supporting R scripts available in skill directory for on-demand reading
 *
 * Usage:
 *   node coding_agent_skills_runner.js --agent claude
 *   node coding_agent_skills_runner.js --agent codex
 *   node coding_agent_skills_runner.js --agent all
 *   node coding_agent_skills_runner.js --agent claude --prompts simple_01,simple_02
 */

import "dotenv/config";
import { exec } from "child_process";
import fs   from "fs";
import path from "path";
import { fileURLToPath } from "url";
import { createHash } from "crypto";

const __dirname = (() => {
  try { return path.dirname(fileURLToPath(import.meta.url)); }
  catch { return process.cwd(); }
})();

const CONFIG = {
  recipes_dir  : path.resolve(__dirname, "..", "recipes"),
  recipe_file  : path.join(__dirname, "..", "recipes", "recipe.md"),
  skill_file   : path.join(__dirname, "..", "recipes", "SKILL.md"),
  dataset_dir  : process.env.BENCH_DATASET_DIR ||
                 "D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19",
  prompts_file : path.join(__dirname, "prompts", "coding_agent_prompts.jsonl"),
  grader_script: path.join(__dirname, "graders", "graders.R"),
  results_dir  : path.join(__dirname, "results"),
  work_base    : path.join(__dirname, "workdirs"),
  r_exec       : process.env.R_EXEC || "C:/Program Files/R/R-4.5.1/bin/Rscript.exe",
  k                     : Number(process.env.BENCH_K || 1),
  timeout_ms            : Number(process.env.BENCH_TIMEOUT_MS   || 600000),
  r_timeout_ms          : Number(process.env.BENCH_R_TIMEOUT_MS || 600000),
  context_window_tokens : Number(process.env.BENCH_CONTEXT_WINDOW_TOKENS || 128000),
  reserved_output_tokens: Number(process.env.BENCH_MAX_TOKENS            || 12000),
  agents: {
    claude: { cmd: () => `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | claude --print"` },
    codex: {
      cmd: () => {
        const model = process.env.CODEX_MODEL || "gpt-5.6-terra";
        return `powershell -NoProfile -Command "Get-Content -Raw final_prompt.txt | codex exec -m ${model} --sandbox read-only --skip-git-repo-check --ephemeral --color never -"`;
      },
    },
  },
};

const run_id = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
function log(msg) { console.log(`[${new Date().toISOString()}] ${msg}`); }

function sha256File(filepath) {
  return createHash("sha256").update(fs.readFileSync(filepath)).digest("hex");
}

function assertSupportedAgent(agent_name) {
  if (!Object.keys(CONFIG.agents).includes(agent_name)) {
    throw new Error(`Skills runner supports only: ${Object.keys(CONFIG.agents).join(", ")}. Received: ${agent_name}`);
  }
}

function cleanPromptWhitespace(text) {
  return String(text || "").replace(/\r\n/g, "\n").replace(/[ \t]+$/gm, "").replace(/\n{3,}/g, "\n\n").trim();
}

function estimateTokens(text) { return Math.ceil(String(text || "").length / 3.5); }

function textStats(text) {
  const v = String(text || "");
  return { chars: v.length, bytes: Buffer.byteLength(v, "utf8"), estimated_tokens: estimateTokens(v), lines: v.length === 0 ? 0 : v.split(/\r?\n/).length };
}

function loadJSONL(filepath) {
  log(`Loading prompts from: ${filepath}`);
  if (!fs.existsSync(filepath)) { log(`ERROR: not found: ${filepath}`); return []; }
  const lines = fs.readFileSync(filepath, "utf8").split("\n").filter(l => l.trim());
  log(`Found ${lines.length} lines in prompts file`);
  const parsed = lines.map(l => { try { return JSON.parse(l); } catch(e) { return null; } }).filter(Boolean);
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
  text = String(text).replace(/<scratchpad>[\s\S]*?<\/scratchpad>/gi, "").replace(/<thinking>[\s\S]*?<\/thinking>/gi, "").trim();
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

function computePassAtK(attempts) { return attempts.some(a => a.grader_pass) ? 1.0 : 0.0; }

function installNativeSkill(agent, work_dir) {
  let skill_dir;
  if (agent === "claude")     skill_dir = path.join(work_dir, ".claude", "skills", "single-cell-r");
  else if (agent === "codex") skill_dir = path.join(work_dir, ".agents", "skills", "single-cell-r");
  else throw new Error(`installNativeSkill: unsupported agent ${agent}`);

  if (!fs.existsSync(CONFIG.skill_file))  throw new Error(`SKILL.md not found: ${CONFIG.skill_file}`);
  if (!fs.existsSync(CONFIG.recipe_file)) throw new Error(`recipe.md not found: ${CONFIG.recipe_file}`);

  fs.mkdirSync(skill_dir, { recursive: true });

  // Copy entire recipes folder (SKILL.md, recipe.md, code/, workflow JSONs)
  // SKILL.md already exists in recipes/ — no overwrite needed
  fs.cpSync(CONFIG.recipes_dir, skill_dir, { recursive: true, force: true });

  // Verify SKILL.md was copied correctly
  const sourceHash    = sha256File(CONFIG.skill_file);
  const installedHash = sha256File(path.join(skill_dir, "SKILL.md"));
  if (sourceHash !== installedHash) throw new Error("Installed SKILL.md does not match source SKILL.md");

  log(`  [skill] installed ${skill_dir}`);
  log(`  [skill] SKILL.md sha256=${installedHash}`);
  return { skill_dir, hash: installedHash };
}

function cleanWorkDir(work_dir) {
  if (!fs.existsSync(work_dir)) return;
  for (const entry of fs.readdirSync(work_dir, { withFileTypes: true }))
    fs.rmSync(path.join(work_dir, entry.name), { recursive: true, force: true });
}

function setupWorkDir(agent, prompt_id, attempt) {
  const work_dir = path.join(CONFIG.work_base, agent, `${prompt_id}_attempt${attempt}`);
  fs.mkdirSync(work_dir, { recursive: true });
  cleanWorkDir(work_dir);
  installNativeSkill(agent, work_dir);
  return work_dir;
}

function buildPromptWithContext(prompt_obj, agent_name) {
  const activation = agent_name === "claude"
    ? `/single-cell-r`
    : `Use the project skill named "single-cell-r" for this task.`;

  const promptText = cleanPromptWhitespace(`${activation}

You are an expert bioinformatician working with single-cell RNA-seq data in R using Seurat.

Use the loaded skill instructions and inspect only the supporting files needed for this specific task.

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
- Generate the smallest complete executable script.`);

  const act_type = agent_name === "claude" ? "slash_command" : "explicit_skill_name";
  log(`  [skill] mode=native_skill_oracle_activated activation=${act_type}`);

  return {
    promptText,
    contextMetadata: {
      mode                        : "native_skill_oracle_activated",
      skill_name                  : "single-cell-r",
      activation                  : act_type,
      delivery                    : "native_progressive_disclosure",
      runner_selected_scripts     : false,
      recipe_injected_into_prompt : false,
      final_prompt                : textStats(promptText),
      configured_context_window_tokens: CONFIG.context_window_tokens,
      reserved_output_tokens          : CONFIG.reserved_output_tokens,
      skill_provenance: {
        canonical_skill_source      : CONFIG.skill_file,
        procedural_recipe_source    : CONFIG.recipe_file,
        installed_filename          : "SKILL.md",
        skill_sha256                : sha256File(CONFIG.skill_file),
        recipe_sha256               : sha256File(CONFIG.recipe_file),
        recipe_injected_into_prompt : false,
        contents_modified           : false,
      },
    },
  };
}

function callAgent(agent_name, prompt_text, work_dir) {
  return new Promise((resolve) => {
    const cmd = CONFIG.agents[agent_name].cmd();
    fs.writeFileSync(path.join(work_dir, "final_prompt.txt"), prompt_text, "utf8");
    log(`  [prompt] chars=${prompt_text.length}`);
    log(`  [${agent_name}] calling agent...`);
    exec(cmd, { timeout: CONFIG.timeout_ms, cwd: work_dir, env: { ...process.env }, maxBuffer: 20 * 1024 * 1024 },
      (error, stdout, stderr) => {
        fs.writeFileSync(path.join(work_dir, "agent_stdout.txt"), stdout || "", "utf8");
        fs.writeFileSync(path.join(work_dir, "agent_stderr.txt"), stderr || "", "utf8");
        if (error) log(`  [${agent_name}] agent error: ${error.message?.slice(0, 180)}`);
        resolve({ stdout: stdout || "", stderr: stderr || "", error: error?.message || null });
      });
  });
}

function executeRCode(r_code, work_dir, prompt_id) {
  return new Promise((resolve) => {
    const r_file = path.join(work_dir, `${prompt_id}_generated.R`);
    fs.writeFileSync(r_file, r_code, "utf8");
    log(`  [R] executing ${path.basename(r_file)}`);
    exec(`"${CONFIG.r_exec}" "${r_file}" 2>&1`, { timeout: CONFIG.r_timeout_ms, cwd: work_dir, maxBuffer: 5 * 1024 * 1024 },
      (error, stdout) => {
        fs.writeFileSync(path.join(work_dir, `${prompt_id}_r_stdout.txt`), stdout || "", "utf8");
        resolve({ success: !error || error.code === 0, stdout: stdout || "", error: error?.message || null, r_file });
      });
  });
}

function runGrader(prompt_id, work_dir) {
  return new Promise((resolve) => {
    const cmd = `"${CONFIG.r_exec}" "${CONFIG.grader_script}" "${prompt_id}" "${work_dir}"`;
    exec(cmd, { timeout: 30000, maxBuffer: 1024 * 1024 }, (error, stdout) => {
      try {
        const lines = (stdout || "").trim().split("\n").filter(l => l.startsWith("{"));
        resolve(JSON.parse(lines.pop() || "{}"));
      } catch {
        resolve({ prompt_id, pass: false, score: 0.0, reason: "Grader parse error: " + (error?.message || "unknown"), grounding_rate: 0.0, grounding_fns: [] });
      }
    });
  });
}

async function evaluatePrompt(agent_name, prompt_obj, attempt = 1) {
  const prompt_id = prompt_obj.id;
  log(`\n  [${agent_name}/${prompt_id}] attempt ${attempt}`);

  const work_dir = setupWorkDir(agent_name, prompt_id, attempt);
  const { promptText: prompt_text, contextMetadata } = buildPromptWithContext(prompt_obj, agent_name);
  fs.writeFileSync(path.join(work_dir, "context_metadata.json"), JSON.stringify(contextMetadata, null, 2), "utf8");

  const t_start      = Date.now();
  const agent_result = await callAgent(agent_name, prompt_text, work_dir);
  const t_agent      = Date.now() - t_start;
  log(`  [${agent_name}/${prompt_id}] agent responded in ${(t_agent / 1000).toFixed(1)}s`);

  const base = {
    agent: agent_name, prompt_id,
    category        : prompt_obj.category || prompt_obj.id?.split("_")[0],
    difficulty      : prompt_obj.difficulty,
    attempt, agent_latency_ms: t_agent, context: contextMetadata, work_dir,
  };

  if (agent_result.error && !agent_result.stdout) {
    return { ...base, total_ms: Date.now() - t_start, r_executed: false, r_success: false,
      grader_pass: false, grader_score: 0, grounding_rate: 0, failure_type: "agent_cli_error", reason: agent_result.error };
  }

  const code = extractRCode(agent_result.stdout);
  if (!code) {
    log(`  [${agent_name}/${prompt_id}] no R code extracted`);
    return { ...base, total_ms: Date.now() - t_start, r_executed: false, r_success: false,
      grader_pass: false, grader_score: 0, grounding_rate: 0, failure_type: "no_code_extracted", reason: "No R code in agent output." };
  }

  const r_result = await executeRCode(code, work_dir, prompt_id);
  const t_total  = Date.now() - t_start;

  if (!r_result.success) {
    log(`  [${agent_name}/${prompt_id}] R ✗ (${(t_total / 1000).toFixed(1)}s)`);
    return { ...base, total_ms: t_total, r_executed: true, r_success: false, r_error: r_result.error,
      grader_pass: false, grader_score: 0, grounding_rate: 0,
      failure_type: r_result.error && /timed out|timeout/i.test(r_result.error) ? "timeout" : "r_runtime_error",
      reason: r_result.error };
  }

  log(`  [${agent_name}/${prompt_id}] R ✓ (${(t_total / 1000).toFixed(1)}s)`);
  const gr = await runGrader(prompt_id, work_dir);
  log(`  [${agent_name}/${prompt_id}] grader: pass=${gr.pass} score=${gr.score?.toFixed(2)} grounding=${gr.grounding_rate?.toFixed(2)}`);
  log(`  [${prompt_id}] Pass@1=${gr.pass ? "1.000" : "0.000"}`);

  return { ...base, total_ms: t_total, r_executed: true, r_success: r_result.success, r_error: r_result.error,
    grader_pass: gr.pass || false, grader_score: gr.score || 0, grounding_rate: gr.grounding_rate || 0,
    grounding_fns: gr.grounding_fns || [], failure_type: gr.pass ? "pass" : "grader_fail", reason: gr.reason || "" };
}

async function runBenchmark(agent_name, prompt_ids = null) {
  assertSupportedAgent(agent_name);
  log(`\n${"═".repeat(70)}`);
  log(`Coding Agent Skills Benchmark  run_id=${run_id}`);
  log(`Agent: ${agent_name}  k=${CONFIG.k}  condition: native_skill_oracle_activated`);
  log(`${"═".repeat(70)}`);

  const all_prompts = loadJSONL(CONFIG.prompts_file);
  const prompts = prompt_ids ? all_prompts.filter(p => prompt_ids.includes(p.id)) : all_prompts;
  log(`Loaded ${prompts.length} prompts`);
  fs.mkdirSync(CONFIG.results_dir, { recursive: true });
  fs.mkdirSync(CONFIG.work_base,   { recursive: true });

  const all_results = [];
  for (const prompt_obj of prompts)
    for (let attempt = 1; attempt <= CONFIG.k; attempt++)
      all_results.push(await evaluatePrompt(agent_name, prompt_obj, attempt));

  const prompt_groups = {};
  for (const r of all_results) { if (!prompt_groups[r.prompt_id]) prompt_groups[r.prompt_id] = []; prompt_groups[r.prompt_id].push(r); }
  const prompt_passes  = Object.values(prompt_groups).map(g => computePassAtK(g));
  const pass_rate      = prompt_passes.reduce((s, v) => s + v, 0) / Math.max(prompt_passes.length, 1);
  const r_success_rate = all_results.filter(r => r.r_success).length    / Math.max(all_results.length, 1);
  const avg_score      = all_results.reduce((s, r) => s + (r.grader_score    || 0), 0) / Math.max(all_results.length, 1);
  const avg_grounding  = all_results.reduce((s, r) => s + (r.grounding_rate  || 0), 0) / Math.max(all_results.length, 1);
  const avg_latency    = all_results.reduce((s, r) => s + (r.agent_latency_ms|| 0), 0) / Math.max(all_results.length, 1);

  log(`\n${"─".repeat(70)}`);
  log(`BENCHMARK SUMMARY — ${agent_name}`);
  log(`${"─".repeat(70)}`);
  log(`Pass@${CONFIG.k} (strict):     ${(pass_rate * 100).toFixed(1)}%`);
  log(`Workflow success (R):  ${(r_success_rate * 100).toFixed(1)}%`);
  log(`Avg grader score:      ${avg_score.toFixed(3)}`);
  log(`Avg grounding rate:    ${(avg_grounding * 100).toFixed(1)}%`);
  log(`Avg agent latency:     ${(avg_latency / 1000).toFixed(1)}s`);

  const cats = {}, diffs = {}, fails = {};
  for (const r of all_results) {
    const c = r.category || "unknown"; if (!cats[c]) cats[c] = { pass: 0, total: 0 }; if (r.grader_pass) cats[c].pass++; cats[c].total++;
    const d = r.difficulty || "unknown"; if (!diffs[d]) diffs[d] = { pass: 0, total: 0 }; if (r.grader_pass) diffs[d].pass++; diffs[d].total++;
    const f = r.failure_type || "unknown"; fails[f] = (fails[f] || 0) + 1;
  }
  log(`\nBy category:`);
  for (const [c, s] of Object.entries(cats)) log(`  ${c.padEnd(20)} ${((s.pass / s.total) * 100).toFixed(1)}% (${s.total} prompts)`);
  log(`\nBy difficulty:`);
  for (const [d, s] of Object.entries(diffs)) log(`  ${d.padEnd(10)} ${((s.pass / s.total) * 100).toFixed(1)}% (${s.total} prompts)`);
  log(`\nBy failure type:`);
  for (const [f, n] of Object.entries(fails)) log(`  ${f.padEnd(20)} ${n}`);

  const out_file = path.join(CONFIG.results_dir, `coding_agent_skills_${agent_name}_${run_id}.json`);
  fs.writeFileSync(out_file, JSON.stringify({
    run_id, agent: agent_name,
    condition       : "native_skill_oracle_activated",
    timestamp       : new Date().toISOString(),
    n_prompts       : prompts.length,
    k               : CONFIG.k,
    metrics         : { pass_at_k: pass_rate, r_success_rate, avg_score, avg_grounding, avg_latency_ms: avg_latency },
    skill_provenance: { canonical_skill_source: CONFIG.skill_file, procedural_recipe_source: CONFIG.recipe_file, installed_filename: "SKILL.md", skill_sha256: sha256File(CONFIG.skill_file), recipe_sha256: sha256File(CONFIG.recipe_file), recipe_injected_into_prompt: false, contents_modified: false },
    results         : all_results,
  }, null, 2));

  log(`\nResults saved → ${out_file}`);
  log(`✅ Done.\n${"═".repeat(70)}\n`);
}

const agent_arg   = process.argv.find((a, i) => process.argv[i - 1] === "--agent") || "claude";
const prompts_arg = process.argv.find((a, i) => process.argv[i - 1] === "--prompts");
const prompt_ids  = prompts_arg ? prompts_arg.split(",").map(s => s.trim()) : null;

if (agent_arg === "all") {
  for (const agent of Object.keys(CONFIG.agents)) await runBenchmark(agent, prompt_ids);
} else {
  await runBenchmark(agent_arg, prompt_ids);
}

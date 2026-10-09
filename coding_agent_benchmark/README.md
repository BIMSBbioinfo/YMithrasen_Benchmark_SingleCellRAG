# Coding-agent benchmark

| File | What it does |
|---|---|
| `coding_agent_runner.js` | Runs Claude Code, Codex CLI or Gemini CLI on the 30 tasks. |
| `graders/graders.R` | One grading function per task. `Rscript graders.R <prompt_id> <work_dir>` prints a JSON verdict. |
| `prompts/coding_agent_prompts.jsonl` | The 30 tasks (same texts as the API benchmark). |

## Conditions

| Condition | In the result files | Flag |
|---|---|---|
| With Skills | `oracle_activated_progressive_disclosure` | default |
| Without Skills | `no_recipe` | `--no-recipe` |

With Skills, the prompt has a short teaser pointing to `./recipes/recipe.md` and
`recipes/` is copied into the agent's working directory. Without Skills, that
folder is not there at all.

## Usage

```bash
node coding_agent_benchmark/coding_agent_runner.js --agent claude
node coding_agent_benchmark/coding_agent_runner.js --agent codex --no-recipe
node coding_agent_benchmark/coding_agent_runner.js --agent gemini --prompts simple_01,hard_04
node coding_agent_benchmark/coding_agent_runner.js --agent claude --resume
```

Settings in `.env`: `BENCH_K` (attempts per prompt, default 1),
`CODEX_MODEL`, `GEMINI_MODEL`, `BENCH_DATASET_DIR`.

If the agent CLI itself fails (session or rate limit), the run stops and that
attempt is not counted. Continue with `--resume`.

## benchmark_version

Each result file stores a hash of the prompts, `graders.R`, everything in
`recipes/` and the dataset path, so runs can be matched to the exact benchmark
content. The files in `results/coding_agents/` have `744037ef0c1d`. The files
here give the same hash when `BENCH_DATASET_DIR` is the original path
(`D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19`).

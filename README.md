# Agentic AI for Single-Cell Omics: Benchmark

A benchmark for AI systems that write R/Seurat code for single-cell RNA-seq
analysis. It contains 30 tasks (10 simple, 10 medium, 10 hard) on the PBMC 3k
dataset. Every generated script is run in R and checked automatically.

**Author:** Yoga Mithra Senthil Kumar (Saarland University / Max Delbrück Center Berlin)
**Supervisors:** Altuna Akalin, Alexander Blume, Andreas Keller

## What is compared

| Systems | Conditions | Code |
|---|---|---|
| 7 API models (GPT-5.5, GPT-5.4 Mini, Gemini 3.5 Flash, Gemini 3.1 Flash-Lite, Claude Opus 4.8, Claude Sonnet 5, Claude Haiku 4.5) | Without Skills / With Skills | `api_benchmark/` |
| 3 coding agents (Claude Code, Codex CLI, Gemini CLI) | Without Skills / With Skills | `coding_agent_benchmark/` |
| Agentic RAG (ROLV), with and without the pre-execution planner | – | `api_benchmark/` |

"With Skills" means the system is told that a relevant skill exists and reads
`recipes/` on demand. API models do this through a `view()` tool, coding agents
with their own file tools.

## Layout

```
api_benchmark/            runners and model shims for the API models and Agentic RAG
coding_agent_benchmark/   runner and grader for the coding agents
recipes/                  the skill (recipe.md, workflows, R scripts)
data/pbmc3k/              PBMC 3k count matrix
results/                  result files, summary tables, figures
analysis/                 scripts for the tables and figures
environment/              R sessionInfo of the machine the runs were made on
```

## Setup

Copy `.env.example` to `.env` and add your API keys. Then use one of the two routes.

**Docker (same environment on Windows, macOS and Linux)**

```bash
docker build -t scrna-benchmark .
docker run --rm -it -v "$PWD/.env:/benchmark/.env:ro" \
  -v "$PWD/api_benchmark/results:/benchmark/api_benchmark/results" \
  -v "$PWD/coding_agent_benchmark/results:/benchmark/coding_agent_benchmark/results" \
  scrna-benchmark
```

The image contains R 4.5.1 with the R packages, Node.js 20 and the `claude`,
`codex` and `gemini` CLIs. The second command opens a shell in `/benchmark`,
where the commands below can be run. Result files are written to the two
`results/` folders on your machine. In PowerShell, write `${PWD}` instead of `$PWD`.

**Without Docker**

Install R 4.5.1 and Node.js 20, then:

```bash
Rscript setup_env.R        # installs the R packages
npm install
```

For the coding-agent benchmark, install the `claude`, `codex` and `gemini`
CLIs and log in. `environment/sessionInfo.txt` lists the package versions of
the R 4.5.1 installation the benchmark was developed with.

Agentic RAG needs access to ROLV and a local copy of the ROLV desktop app
(`ROLV_ORCHESTRATOR_PATH` in `.env`) and is not part of the Docker image.

## Running

Run from the repository root.

```bash
# API models, Without Skills (pick any of the model keys)
node api_benchmark/benchmark_runner.js --k 3 --model gpt-5.5,claude-sonnet-5

# API models, With Skills
node api_benchmark/skills_benchmark_runner.js --k 3 --model gpt-5.5,claude-sonnet-5

# Agentic RAG
node api_benchmark/benchmark_runner.js --k 3 --model rolv,rolv-nopreexec,rolv-raw

# Coding agents (BENCH_K in .env sets the attempts per prompt)
node coding_agent_benchmark/coding_agent_runner.js --agent claude
node coding_agent_benchmark/coding_agent_runner.js --agent claude --no-recipe

# Tables and figures
python analysis/summarise_results.py
python analysis/make_figures.py
```

New runs are written to `api_benchmark/results/` and
`coding_agent_benchmark/results/`, which are not tracked by git.

## Metrics

- **Pass@k**: share of prompts where at least one of k attempts passed.
- **Attempt success rate**: share of single attempts that passed.
- **Workflow score**: 0.7 × expected pipeline steps present + 0.3 × steps in the right order.
- **Grounding rate**: share of function calls that exist in the installed R packages or are defined in the script.
- **R success rate / completion rate / grader score** (coding agents): R code ran / attempt passed the grader / mean grader score.

For API models an attempt passes when the code runs and the prompt's validation
rules hold. For coding agents it passes when `graders/graders.R` accepts the output.

## Results

Pass@3 on the 30 prompts, in %. All metrics and the figures are in `results/`.

| API model | Without Skills | With Skills |
|---|---|---|
| GPT-5.5 | 70.0 | 93.3 |
| GPT-5.4 Mini | 70.0 | 90.0 |
| Gemini 3.5 Flash | 86.7 | 96.7 |
| Gemini 3.1 Flash-Lite | 73.3 | 83.3 |
| Claude Opus 4.8 | 76.7 | 90.0 |
| Claude Sonnet 5 | 90.0 | 90.0 |
| Claude Haiku 4.5 | 40.0 | 86.7 |

| Coding agent | Without Skills | With Skills |
|---|---|---|
| Claude Code | 93.3 | 96.7 |
| Codex CLI | 90.0 | 100.0 |
| Gemini CLI | 93.3 | 96.7 |

| Agentic RAG | Pass@3 |
|---|---|
| Agentic RAG | 93.3 |
| Agentic RAG (no pre-exec) | 93.3 |
| Claude Sonnet 4.5 (raw model) | 16.7 |

The Agentic RAG run was graded on execution only, see `results/README.md`.

## Using other models, agents or data

**Run a subset.** `--model`, `--prompts` and `--agent` select which of the
existing models, prompts and agents are run. `CODEX_MODEL` and `GEMINI_MODEL`
in `.env` change the model behind those two CLIs.

**Add an API model.** Models are defined in code, not in a config file:

1. `api_benchmark/scripts/baseline_shim.js`: add a call function and a branch
   for the new key in the model dispatch at the bottom.
2. `api_benchmark/benchmark_runner.js`: add the key to `MODEL_REQUIRED_KEYS`,
   `ALL_SMOKE_TESTS`, `API_MODELS`, `MODEL_DISPLAY_NAMES` and `MODEL_SAMPLING`.
3. For the With Skills condition, do the same in `scripts/skills_shim.js` and
   in `API_MODELS`, `MODEL_DISPLAY_NAMES` and `MODEL_REQUIRED_KEYS` of
   `skills_benchmark_runner.js`.

For a model with an OpenAI-compatible endpoint, copy one of the existing call
functions and change the model name and base URL.

**Add a coding agent.** Add an entry with the CLI command to `CONFIG.agents`
in `coding_agent_runner.js` (and to the list used by `--agent all` at the bottom
of the file).

**Use another dataset.** `BENCH_DATASET_DIR` points the runners at another
folder of 10x files (`barcodes.tsv`, `genes.tsv`, `matrix.mtx`). The
tasks themselves are written for PBMC 3k, though: the prompt texts describe that
dataset, and the validation rules in `benchmark_prompts.jsonl` and the checks in
`graders.R` (cell counts, marker genes, cluster numbers) assume it. A different
dataset therefore needs its own prompts, validation rules and grader functions.

**Change the skill.** Edit the files in `recipes/`. The teaser is built from
the `name` field and the "What this recipe does" section of `recipe.md`.
# Benchmark: Single-Cell RAG and Agentic Skills Evaluation

Master thesis benchmark — Saarland University / MDC Berlin  
**Author:** Yoga Mithra  
**Supervisors:** Alexander Blume, Altuna Akalin

---

## Overview

This repository contains the benchmarking framework for evaluating AI-assisted
single-cell RNA-seq analysis. It measures two research questions:

1. **RAG vs Vanilla API** — Does a RAG-augmented system (ROLV) generate better
   scRNA-seq code than a vanilla LLM API with no domain knowledge?

2. **Augmented CLI vs Vanilla CLI** — Does providing a coding agent CLI with a
   structured recipe/skill file improve code generation compared to the agent
   working from scratch?

---

## Repository structure

```
benchmark/
  benchmark_runner.js        Main benchmark harness
  prompts/
    benchmark_prompts.jsonl  7 PBMC code generation prompts
    routing_queries.jsonl    30 labelled routing queries
    genomics_qa.jsonl        20 genomics QA pairs
  scripts/
    orchestrator_shim.js     Bridges to ROLV orchestrator.js (API suite)
    baseline_shim.js         GPT/Gemini/Groq API callers
    cli_runner.js            Runs Gemini/Codex/Claude Code CLIs (CLI suite)
    bio_eval.R               Biological metrics (ARI, kBET, LISI, ASW)
    plot_results.R           Publication-ready thesis figures
  datasets/
    filtered_gene_bc_matrices/hg19/   PBMC 3k 10x data

recipes/
  scrna-seurat-pbmc/
    recipe.md                Human-readable workflow description
    workflow.json            Machine-readable workflow graph
    env/r-packages.txt       R dependencies
    code/                    R scripts for each pipeline step
```

---

## Requirements

### Node.js
- Node.js >= 18
- Install dependencies: `cd benchmark && npm install`

### R
- R >= 4.1.0
- Seurat >= 5.0.0
- Install packages listed in `recipes/scrna-seurat-pbmc/env/r-packages.txt`

### API keys (for API suite)
Add to `benchmark/.env`:
```
ROLV_API_URL=https://api.rolv.io/api/chat/completions
ROLV_JWT_TOKEN=<your token>
GROQ_API_KEY=<your key — free at console.groq.com/keys>
GOOGLE_API_KEY=<your key>
R_EXEC=C:/Program Files/R/R-4.5.1/bin/Rscript.exe
```

### CLI tools (for CLI suite)
Install at least one:
- Gemini CLI: `npm install -g @google/gemini-cli`
- Claude Code: `npm install -g @anthropic-ai/claude-code`
- OpenAI Codex CLI: `npm install -g @openai/codex`

---

## Running the benchmark

### API suite (RAG vs Vanilla API)
```bash
cd benchmark
node benchmark_runner.js --suite internal        # Pass@k, workflow, grounding
node benchmark_runner.js --suite routing         # Routing accuracy
node benchmark_runner.js --suite qa              # Gene-Turing Score
```

### CLI suite (Augmented vs Vanilla CLI)
```bash
node benchmark_runner.js --suite cli             # All CLI agents with and without recipe
node benchmark_runner.js --suite cli --agent gemini    # Gemini CLI only
node benchmark_runner.js --suite cli --agent claude    # Claude Code only
node benchmark_runner.js --suite cli --agent codex     # Codex CLI only
```

### Generate figures
```bash
Rscript scripts/plot_results.R \
  --csv results/results_latest_summary.csv \
  --out_dir results/figures
```

---

## Results

Results are saved to `benchmark/results/` as:
- `results_<timestamp>.json` — full raw results
- `results_<timestamp>_summary.csv` — summary table for R

---

## Prompts used

All 7 PBMC prompts follow a naturalistic biologist style — each starts from raw
10x Genomics files (barcodes.tsv, genes.tsv, matrix.mtx) and describes the task
in biological terms rather than pre-specifying function names. Difficulty ranges
from simple (load data, QC, clustering) to moderate (cell type annotation,
differential expression with volcano plot).

See `benchmark/prompts/benchmark_prompts.jsonl` for the full prompt set.

---

## Current results (ROLV, May 2026)

| Metric | Score |
|--------|-------|
| Pass@3 | 100% (20/21 attempts) |
| Workflow Success Rate | 77.1% |
| Grounding Rate | 85.5% |
| Gene-Turing Score | 0.864 |
| Routing Accuracy | 86.7% |
| Routing Macro-F1 | 0.860 |

Baseline comparison (Groq, Gemini, ROLV raw) — pending API key setup.
CLI suite (Gemini CLI, Claude Code, Codex) — pending CLI installation.

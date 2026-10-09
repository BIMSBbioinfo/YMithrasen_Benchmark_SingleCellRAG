# API benchmark

| File | What it does |
|---|---|
| `benchmark_runner.js` | Without Skills baseline and Agentic RAG. |
| `skills_benchmark_runner.js` | With Skills: skill teaser plus a `view()` tool over `recipes/`. |
| `scripts/baseline_shim.js` | One model call, returns the extracted R code. |
| `scripts/skills_shim.js` | Same, with the `view()` tool loop. |
| `scripts/orchestrator_shim.js` | Calls the ROLV orchestrator (`rolv`). |
| `scripts/orchestrator_shim_nopreexec.js` | Same without the planner (`rolv-nopreexec`). |
| `prompts/benchmark_prompts.jsonl` | 30 tasks with `validation` and `expected_workflow`. |

## One attempt

1. The runner starts a shim, which calls the model and returns R code.
2. The code is run with `Rscript --vanilla` inside the dataset folder (10 min limit),
   followed by the prompt's validation checks.
3. The attempt passes if the code ran and all checks held.
4. Workflow score and grounding rate are computed from the same code.

Before a run, the runner checks keys, R packages, dataset and prompts and makes
one test call per model.

## Model keys

| Key | Model |
|---|---|
| `gpt-5.5`, `gpt-5.4-mini` | OpenAI |
| `gemini-openrouter`, `gemini-flash-lite` | Gemini 3.5 Flash, Gemini 3.1 Flash-Lite via OpenRouter |
| `claude-opus-4-8`, `claude-sonnet-5`, `claude-haiku-4-5` | Anthropic |
| `rolv`, `rolv-nopreexec`, `rolv-raw` | Agentic RAG, without planner, raw model (only in `benchmark_runner.js`) |

GPT-5.4 Mini and the Gemini models run at temperature 0.6. GPT-5.5 and the
Claude models do not take a temperature.

## Options

`benchmark_runner.js`: `--k`, `--model`, `--prompts`, `--allow-partial-prompts`

`skills_benchmark_runner.js`: `--k`, `--model`, `--prompts`, `--allow-partial-prompts`,
`--baseline <file>`, `--skip-comparison`

With Skills, a model can explore the skill folder for up to 12 turns
(`BENCH_MAX_TOOL_TURNS`) and then gets one last turn without tools.

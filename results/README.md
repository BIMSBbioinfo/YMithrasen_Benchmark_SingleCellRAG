# Results

All runs use the 30 prompts. k = 3 attempts per prompt unless stated.

## Files

**API models and Agentic RAG (`api/`)**

| File | What |
|---|---|
| `vanilla_k3_merged.json` | Without Skills, 7 models |
| `skills_results_2026-08-02T13-16-53-819.json` | With Skills, 7 models |
| `results_2026-06-11T20-49-17-589.json` | Agentic RAG, Agentic RAG without pre-exec, Claude Sonnet 4.5 (raw model) |
| `results_2026-07-19T20-16-10-541.json` | Without Skills, earlier k = 3 run (superseded by the merged file) |
| `results_2026-07-19T16-55-57-401.json` | Without Skills, k = 1 |

**Coding agents (`coding_agents/`)**

| File | What |
|---|---|
| `coding_agent_<agent>_<run>.json` | Without Skills (Claude Code, Codex CLI, Gemini CLI) |
| `coding_agent_skills-<agent>_<run>.json` | With Skills |
| `token_usage/` | usage exports for the Gemini CLI and Codex runs (account ids removed) |

**Tables and figures**

- `api_summary.csv`, `coding_agents_summary.csv`: one row per model or agent
  and condition. Written by `analysis/summarise_results.py`.
- `figures/`: the figures, and `stats_table.csv` with the mean, SD and 95%
  interval behind every bar (Wilson for Pass@3, task-level bootstrap otherwise).
  `analysis/make_figures.py` redraws the figures.

## Reading the result files

API files:

```
results.<model>.passK.avg                     Pass@k, mean over prompts
results.<model>.passK.avgAttemptSuccessRate   share of attempts that passed
results.<model>.passK.results[]               per prompt: score and every attempt (code, passed, error)
results.<model>.quality.avgWorkflow           workflow score
results.<model>.quality.avgGrounding          grounding rate
```

The Skills file also stores the `view()` calls of each attempt
(`tool_calls_made`) and whether the turn limit was reached
(`tool_budget_exhausted`).

Coding-agent files:

```
metrics.pass_at_k        share of prompts with at least one passing attempt
metrics.r_success_rate   share of attempts whose R code ran
metrics.avg_score        mean grader score
metrics.avg_latency_ms   mean agent response time
by_difficulty            Pass@k for simple, medium and hard prompts
results[]                one record per attempt (grader_pass, grader_score, failure_type, reason)
```

## Main findings

- **Skills help the API models.** Pass@3 rises for six of seven models and
  stays the same for Claude Sonnet 5. The largest gain is Claude Haiku 4.5
  (40.0% to 86.7%). Workflow and grounding are already above 90% without
  Skills and change little.
- **The gain is largest on hard tasks.** Mean attempt success across the seven
  models goes from 81.4% to 99.0% on simple, 66.7% to 91.0% on medium and
  38.1% to 58.1% on hard prompts.
- **Coding agents start high and still improve.** Pass@3 is 90–93% without
  Skills and 97–100% with Skills. Codex CLI gains most (completion rate 74.4%
  to 94.4%).
- **Skills cost time and tokens for agents.** Mean latency goes up for all
  three agents, and total tokens rise by 1.85× (Claude Code), 5.54× (Codex CLI)
  and 8.46× (Gemini CLI).
- **Agentic RAG** reaches 93.3% Pass@3 with and without the pre-execution
  planner, against 16.7% for the same model without retrieval.

## Notes on the files

- `vanilla_k3_merged.json` combines two runs: Claude Haiku 4.5 comes from a
  re-run, because 16 of its 90 attempts in the first run failed with connection
  errors. The file's `merge_provenance` field has the details.
- Groq (GPT-OSS 20B) was dropped from the benchmark, and its entry was removed
  from the files that contained it. Those files have a `post_processing` field
  with the checksum of the original.
- The Agentic RAG file was graded on execution only: an attempt passed if the R
  code ran without error, with no task validation. Its numbers are therefore
  not on the same scale as the Without/With Skills files. The `grading` column
  of `api_summary.csv` shows which grading each file used.
- The API Skills run recorded a different `recipe.md` (`skill_provenance.raw_sha256`
  starts with `07309b3a`) than the one in this repository, which is the version
  the coding-agent runs used.
- Paths inside the files are from the machine the runs were made on.

"""
make_figures.py

Draws the figures in results/figures/ from the result files.

Bar heights are read from the JSON files in results/. Error bars are the
95% intervals stored in results/figures/stats_table.csv (Wilson for Pass@3,
task-level bootstrap for everything else).

Usage (from the repository root):
    python analysis/make_figures.py [output_dir]

Figures are written as 300-dpi PNG to output_dir (default: results/figures/).
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path(__file__).resolve().parent.parent / "results"
STATS_TABLE = RESULTS / "figures" / "stats_table.csv"
OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS / "figures"

VANILLA_FILE = RESULTS / "api" / "vanilla_k3_merged.json"
SKILLS_FILE = RESULTS / "api" / "skills_results_2026-08-02T13-16-53-819.json"
RAG_FILE = RESULTS / "api" / "results_2026-06-11T20-49-17-589.json"
CODING_FILES = {
    ("claude", "without"): "coding_agent_claude_2026-08-16T14-23-33.json",
    ("claude", "with"): "coding_agent_skills-claude_2026-08-16T11-20-54.json",
    ("codex", "without"): "coding_agent_codex_2026-08-16T20-33-59.json",
    ("codex", "with"): "coding_agent_skills-codex_2026-08-16T17-09-28.json",
    ("gemini", "without"): "coding_agent_gemini_2026-08-17T12-23-36.json",
    ("gemini", "with"): "coding_agent_skills-gemini_2026-08-17T08-54-26.json",
}

# model order used on the x-axis of the API figures
MODELS = ["claude-sonnet-5", "gemini-flash-lite", "gpt-5.4-mini", "gemini-openrouter",
          "gpt-5.5", "claude-opus-4-8", "claude-haiku-4-5"]
MODEL_NAMES = {
    "gpt-5.5": "GPT-5.5",
    "gpt-5.4-mini": "GPT-5.4 Mini",
    "gemini-openrouter": "Gemini 3.5 Flash",
    "gemini-flash-lite": "Gemini 3.1 Flash-Lite",
    "claude-opus-4-8": "Claude Opus 4.8",
    "claude-sonnet-5": "Claude Sonnet 5",
    "claude-haiku-4-5": "Claude Haiku 4.5",
}
AGENTS = ["claude", "codex", "gemini"]
AGENT_NAMES = {"claude": "Claude Code", "codex": "Codex", "gemini": "Gemini CLI"}
RAG_SYSTEMS = ["rolv", "rolv-nopreexec", "rolv-raw"]
RAG_NAMES = {"rolv": "Agentic RAG", "rolv-nopreexec": "Agentic RAG without PreExec",
             "rolv-raw": "Claude Sonnet 4.5 (Raw)"}
TIERS = ["simple", "medium", "hard"]
CONDITIONS = ["without", "with"]
CONDITION_NAMES = {"without": "Without Skills", "with": "With Skills"}

GREY, RED = "#8c8c8c", "#c0392b"
CONDITION_COLOURS = {"without": GREY, "with": RED}
NAVY, BLUE, ORANGE, GREEN = "#1b3a5c", "#5b9bd5", "#eba63a", "#8fae53"

# total tokens per run, from the provider usage records
TOKENS = {
    "claude": {"without": 46245, "with": 85579},
    "codex": {"without": 14566, "with": 80703},
    "gemini": {"without": 23719, "with": 200701},
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
})


# ---- loading ----------------------------------------------------------------

def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def load_intervals():
    """stats_table.csv as {key: (ci_low, ci_high)}."""
    with open(STATS_TABLE, newline="", encoding="utf-8") as fh:
        return {r["key"]: (float(r["ci_low"]), float(r["ci_high"])) for r in csv.DictReader(fh)}


def api_metrics(data):
    """Summary metrics in percent, per model."""
    out = {}
    for model in MODELS:
        entry = data["results"][model]
        out[model] = {
            "pass3": entry["passK"]["avg"] * 100,
            "asr": entry["passK"]["avgAttemptSuccessRate"] * 100,
            "workflow": entry["quality"]["avgWorkflow"] * 100,
            "grounding": entry["quality"]["avgGrounding"] * 100,
        }
    return out


def api_by_tier(data):
    """Mean attempt success rate (percent) per model and difficulty tier."""
    out = {}
    for model in MODELS:
        prompts = data["results"][model]["passK"]["results"]
        for tier in TIERS:
            vals = [p["attempt_success_rate"] for p in prompts if p["id"].startswith(tier)]
            out[(model, tier)] = np.mean(vals) * 100
    return out


def coding_metrics():
    out = {}
    for (agent, cond), name in CODING_FILES.items():
        data = load_json(RESULTS / "coding_agents" / name)
        attempts = data["results"]
        out[(agent, cond)] = {
            "pass3": data["metrics"]["pass_at_k"] * 100,
            "r_success": data["metrics"]["r_success_rate"] * 100,
            "grader_score": data["metrics"]["avg_score"] * 100,
            "completion": np.mean([bool(a["grader_pass"]) for a in attempts]) * 100,
            "latency": data["metrics"]["avg_latency_ms"] / 1000,
            "by_tier": {t: data["by_difficulty"][t] * 100 for t in TIERS},
        }
    return out


# ---- drawing helpers --------------------------------------------------------

def paired_bars(ax, values, intervals=None, fmt="{:.1f}", fontsize=10, width=0.35, label_gap=1.5):
    """Grey/red bar pair per group. values[cond] is a list, one entry per group."""
    x = np.arange(len(values["without"]))
    for i, cond in enumerate(CONDITIONS):
        pos = x + (i - 0.5) * (width + 0.01)
        vals = np.array(values[cond], dtype=float)
        ax.bar(pos, vals, width=width, color=CONDITION_COLOURS[cond], label=CONDITION_NAMES[cond])
        tops = vals
        if intervals is not None:
            lo = np.array([c[0] for c in intervals[cond]])
            hi = np.array([c[1] for c in intervals[cond]])
            ax.errorbar(pos, vals, yerr=[np.clip(vals - lo, 0, None), np.clip(hi - vals, 0, None)],
                        fmt="none", ecolor="black", elinewidth=1, capsize=3)
            tops = hi
        for p, v, top in zip(pos, vals, tops):
            ax.text(p, top + label_gap, fmt.format(v), ha="center", va="bottom", fontsize=fontsize)
    ax.set_xticks(x)
    return x


def metric_bars(ax, groups, metrics, values, intervals=None, fmt="{:.0f}", width=0.26, fontsize=10):
    """Several metrics side by side per group. metrics is a list of (key, label, colour)."""
    x = np.arange(len(groups))
    for i, (key, label, colour) in enumerate(metrics):
        pos = x + (i - (len(metrics) - 1) / 2) * (width + 0.01)
        vals = np.array([values[g][key] for g in groups], dtype=float)
        ax.bar(pos, vals, width=width, color=colour, label=label)
        tops = vals
        if intervals is not None:
            lo = np.array([intervals[(g, key)][0] for g in groups])
            hi = np.array([intervals[(g, key)][1] for g in groups])
            ax.errorbar(pos, vals, yerr=[np.clip(vals - lo, 0, None), np.clip(hi - vals, 0, None)],
                        fmt="none", ecolor="black", elinewidth=1, capsize=3)
            tops = hi
        for p, v, top in zip(pos, vals, tops):
            ax.text(p, top + 1.5, fmt.format(v), ha="center", va="bottom", fontsize=fontsize)
    ax.set_xticks(x)


def finish(fig, ax, title, ncol=2, top=0.86, legend_y=0.94, title_y=1.11):
    """Fit the axes, then put the legend and the title above them."""
    for a in fig.axes:
        a.tick_params(labelsize=12)
        a.yaxis.label.set_size(12)
    fig.tight_layout(rect=[0, 0, 1, top])
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=ncol, bbox_to_anchor=(0.5, legend_y),
               frameon=False, fontsize=12)
    fig.suptitle(title, y=title_y, fontsize=14)


def save(fig, name):
    fig.savefig(OUTDIR / f"{name}.png", bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


# ---- API models -------------------------------------------------------------

def fig_api_metric(api, ci, metric, ylabel, title, name):
    """One metric per model, with and without Skills, with 95% intervals."""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    values = {c: [api[c][m][metric] for m in MODELS] for c in CONDITIONS}
    intervals = {c: [ci[f"api|{m}|{c}|{metric}"] for m in MODELS] for c in CONDITIONS}
    paired_bars(ax, values, intervals)
    ax.set_xticklabels([MODEL_NAMES[m] for m in MODELS], rotation=20, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 112)
    finish(fig, ax, title)
    save(fig, name)


def fig_api_grounding(api):
    """Grounding rate per model, point estimates only."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(MODELS))
    for i, (cond, colour) in enumerate([("without", "#888888"), ("with", "#d62728")]):
        vals = [api[cond][m]["grounding"] for m in MODELS]
        pos = x + (i - 0.5) * 0.35
        ax.bar(pos, vals, width=0.35, color=colour, edgecolor="white", label=CONDITION_NAMES[cond])
        for p, v in zip(pos, vals):
            ax.text(p, v + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_NAMES[m] for m in MODELS], rotation=20, ha="right")
    ax.set_ylabel("Grounding rate (%)")
    ax.set_ylim(0, 112)
    ax.set_title("API models: Without vs. With Skills — Grounding rate", pad=45)
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.14), frameon=False)
    fig.tight_layout()
    save(fig, "fig_api_vanilla_vs_skills_grounding")


def fig_skills_k3(api):
    """Pass@3, workflow and grounding with Skills, models sorted by Pass@3."""
    in_file_order = list(load_json(SKILLS_FILE)["results"])
    order = sorted(in_file_order, key=lambda m: -api["with"][m]["pass3"])
    metrics = [("pass3", "Pass@3", "#1f3864"), ("workflow", "Workflow score", "#eda83c"),
               ("grounding", "Grounding rate", GREEN)]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(order))
    for i, (key, label, colour) in enumerate(metrics):
        vals = [api["with"][m][key] for m in order]
        pos = x + (i - 1) * 0.27
        ax.bar(pos, vals, width=0.27, color=colour, edgecolor="white", linewidth=0.5, label=label)
        for p, v in zip(pos, vals):
            ax.text(p, v + 1.5, f"{v:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_NAMES[m] for m in order], rotation=15, ha="right", fontsize=10.5)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 112)
    ax.set_title("Recipe injection benchmark — Pass@3, Workflow score, Grounding rate", pad=38, fontsize=12)
    ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.13), frameon=False, fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig_skills_k3")


def fig_api_by_difficulty(tiers, ci):
    """Attempt success rate per difficulty tier, averaged over the seven models."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    values = {c: [np.mean([tiers[c][(m, t)] for m in MODELS]) for t in TIERS] for c in CONDITIONS}
    intervals = {c: [ci[f"api_by_diff|{t}|{c}"] for t in TIERS] for c in CONDITIONS}
    paired_bars(ax, values, intervals, width=0.34, label_gap=2)
    ax.set_xticklabels([t.capitalize() for t in TIERS])
    ax.set_ylabel("Attempt success rate (%)")
    ax.set_ylim(0, 115)
    finish(fig, ax, "API-model attempt success by task difficulty\n(mean across 7 models; 95% task-level "
           "bootstrap CI; 10 tasks per difficulty, 7 fixed models)", top=0.88)
    save(fig, "fig_api_by_difficulty")


def fig_api_by_difficulty_per_model(tiers, ci):
    """Same as above, one panel per model (four on top, three centred below)."""
    fig = plt.figure(figsize=(20, 8.5))
    grid = fig.add_gridspec(2, 8, left=0.04, right=0.99, bottom=0.05, top=0.86, hspace=0.55, wspace=0.35)
    slots = [grid[0, 0:2], grid[0, 2:4], grid[0, 4:6], grid[0, 6:8],
             grid[1, 1:3], grid[1, 3:5], grid[1, 5:7]]
    first = None
    for model, slot in zip(MODELS, slots):
        ax = fig.add_subplot(slot)
        first = first or ax
        x = np.arange(len(TIERS))
        for i, cond in enumerate(CONDITIONS):
            vals = np.array([tiers[cond][(model, t)] for t in TIERS])
            lo, hi = zip(*[ci[f"api_by_diff_model|{model}|{t}|{cond}"] for t in TIERS])
            pos = x + (i - 0.5) * 0.36
            ax.bar(pos, vals, width=0.34, color=CONDITION_COLOURS[cond], label=CONDITION_NAMES[cond])
            ax.errorbar(pos, vals, yerr=[np.clip(vals - lo, 0, None), np.clip(np.array(hi) - vals, 0, None)],
                        fmt="none", ecolor="black", elinewidth=0.8, capsize=2)
            for p, v in zip(pos, vals):
                ax.text(p, v + 2, f"{v:.0f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels([t.capitalize() for t in TIERS])
        ax.set_ylim(0, 118)
        ax.set_yticks(range(0, 101, 20))
        ax.set_title(MODEL_NAMES[model], fontsize=12)
        ax.tick_params(labelsize=10.5)
        if model in (MODELS[0], MODELS[4]):
            ax.set_ylabel("Attempt success rate (%)", fontsize=11)
    handles, labels = first.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.95), frameon=False, fontsize=11)
    fig.suptitle("API-model attempt success by task difficulty, with and without procedural Skills "
                 "(per model)\n95% bootstrap CI (n=10 tasks per bar)", y=1.055, fontsize=13)
    save(fig, "fig_api_by_difficulty_per_model")


# ---- Agentic RAG ------------------------------------------------------------

def fig_rag(ci):
    """Agentic RAG, with and without the planner, against the raw model."""
    data = load_json(RAG_FILE)["results"]
    values = {s: {"pass3": data[s]["passK"]["avg"] * 100,
                  "workflow": data[s]["quality"]["avgWorkflow"] * 100,
                  "grounding": data[s]["quality"]["avgGrounding"] * 100} for s in RAG_SYSTEMS}
    intervals = {(s, k): ci[f"rag|{s}|{k}"] for s in RAG_SYSTEMS for k in ("pass3", "workflow", "grounding")}
    metrics = [("pass3", "Pass@3", NAVY), ("workflow", "Workflow score", BLUE),
               ("grounding", "Grounding rate", ORANGE)]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    metric_bars(ax, RAG_SYSTEMS, metrics, values, intervals, fmt="{:.1f}", width=0.245)
    ax.set_xticklabels([RAG_NAMES[s] for s in RAG_SYSTEMS])
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 112)
    finish(fig, ax, "Agentic RAG vs. Raw-Model Baseline\nPass@3 (95% Wilson CI); Workflow score, "
           "Grounding rate (95% task-level bootstrap CI)", ncol=3, top=0.84, legend_y=0.925)
    save(fig, "fig_rag_vs_vanilla_k3")


# ---- coding agents ----------------------------------------------------------

def fig_coding_eval(coding):
    """Pass@3, R success rate and grader score per agent and condition."""
    groups = [(a, c) for a in AGENTS for c in CONDITIONS]
    metrics = [("pass3", "Pass@3", NAVY), ("r_success", "R execution success rate", BLUE),
               ("grader_score", "Deterministic grader score", ORANGE)]
    fig, ax = plt.subplots(figsize=(12.6, 5.8))
    metric_bars(ax, groups, metrics, coding, fmt="{:.0f}", width=0.25)
    ax.set_xticklabels([f"{AGENT_NAMES[a]}\n({CONDITION_NAMES[c]})" for a, c in groups])
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 112)
    finish(fig, ax, "Coding-agent performance with and without procedural Skills\n"
           "(point estimates)", ncol=3, top=0.85, title_y=1.10)
    save(fig, "fig_coding_agent_eval_metrics")


def fig_coding_by_difficulty(coding, ci):
    """Pass@3 per difficulty tier, one panel per agent."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, agent in zip(axes, AGENTS):
        values = {c: [coding[(agent, c)]["by_tier"][t] for t in TIERS] for c in CONDITIONS}
        intervals = {c: [ci[f"coding_by_diff|{agent}|{t}|{c}"] for t in TIERS] for c in CONDITIONS}
        paired_bars(ax, values, intervals, fmt="{:.0f}", width=0.34)
        ax.set_xticklabels([t.capitalize() for t in TIERS])
        ax.set_title(AGENT_NAMES[agent], fontsize=14)
        ax.set_ylim(0, 118)
    axes[0].set_ylabel("Pass rate (%)")
    finish(fig, axes[0], "Coding-agent Pass@3 by task difficulty, with and without procedural Skills\n"
           "95% Wilson CI (n=10 tasks per bar)", top=0.85, legend_y=0.95, title_y=1.13)
    save(fig, "fig_coding_agent_by_difficulty")


def fig_coding_combined(coding, ci):
    """Completion rate, latency and total tokens per agent."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    names = [AGENT_NAMES[a] for a in AGENTS]

    panels = [("completion", "completion", "Completion rate (%)", "Completion rate", (0, 112), 1.5),
              ("latency", "latency_s", "Latency (s)", "Mean agent-response latency", None, 1.2)]
    for ax, (key, ci_key, ylabel, title, ylim, gap) in zip(axes, panels):
        values = {c: [coding[(a, c)][key] for a in AGENTS] for c in CONDITIONS}
        intervals = {c: [ci[f"coding|{a}|{c}|{ci_key}"] for a in AGENTS] for c in CONDITIONS}
        paired_bars(ax, values, intervals, width=0.34, label_gap=gap)
        ax.set_xticklabels(names, rotation=15, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=14)
        if ylim:
            ax.set_ylim(*ylim)

    ax = axes[2]
    x = np.arange(len(AGENTS))
    highest = max(TOKENS[a]["with"] for a in AGENTS)
    for i, cond in enumerate(CONDITIONS):
        vals = [TOKENS[a][cond] for a in AGENTS]
        pos = x + (i - 0.5) * 0.35
        ax.bar(pos, vals, width=0.34, color=CONDITION_COLOURS[cond])
        for p, v in zip(pos, vals):
            ax.text(p, v + highest * 0.01, f"{v:,}", ha="center", va="bottom", fontsize=10)
    for i, agent in enumerate(AGENTS):
        ratio = TOKENS[agent]["with"] / TOKENS[agent]["without"]
        ax.text(i, TOKENS[agent]["with"] * 1.12, f"{ratio:.2f}×", ha="center", va="bottom",
                fontsize=11, fontweight="bold")
    ax.set_ylim(0, highest * 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("Total tokens")
    ax.set_title("Total tokens (descriptive; no per-task\ntoken record in these runs)", fontsize=14)

    finish(fig, axes[0], "Coding-agent reliability and computational overhead with and without procedural "
           "Skills\n(95% bootstrap CI on completion rate and latency; token totals are descriptive run totals)",
           top=0.80, title_y=1.10)
    save(fig, "fig_coding_agent_combined_metrics")


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ci = load_intervals()
    vanilla, skills = load_json(VANILLA_FILE), load_json(SKILLS_FILE)
    api = {"without": api_metrics(vanilla), "with": api_metrics(skills)}
    tiers = {"without": api_by_tier(vanilla), "with": api_by_tier(skills)}
    coding = coding_metrics()

    fig_api_metric(api, ci, "pass3", "Pass@3 (%)",
                   "API models: Pass@3 with and without procedural Skills\n(95% Wilson CI)",
                   "fig_api_vanilla_vs_skills_pass_k")
    fig_api_metric(api, ci, "asr", "Attempt success rate (%)",
                   "API models: attempt success rate with and without procedural Skills\n"
                   "(95% bootstrap CI, task-level)",
                   "fig_api_vanilla_vs_skills_attempt_success_rate")
    fig_api_metric(api, ci, "workflow", "Workflow score (%)",
                   "API models: workflow score with and without procedural Skills\n"
                   "(95% bootstrap CI, task-level)",
                   "fig_api_vanilla_vs_skills_workflow")
    fig_api_grounding(api)
    fig_skills_k3(api)
    fig_api_by_difficulty(tiers, ci)
    fig_api_by_difficulty_per_model(tiers, ci)
    fig_rag(ci)
    fig_coding_eval(coding)
    fig_coding_by_difficulty(coding, ci)
    fig_coding_combined(coding, ci)

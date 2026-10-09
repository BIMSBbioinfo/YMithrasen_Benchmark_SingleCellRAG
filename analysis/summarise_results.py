"""
summarise_results.py

Writes results/api_summary.csv and results/coding_agents_summary.csv from the
result files in results/.

Usage (from the repository root):
    python analysis/summarise_results.py
"""

import csv
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

# (file in results/api/, condition)
API_FILES = [
    ("vanilla_k3_merged.json", "Without Skills"),
    ("skills_results_2026-08-02T13-16-53-819.json", "With Skills"),
    ("results_2026-06-11T20-49-17-589.json", "Agentic RAG"),
    ("results_2026-07-19T20-16-10-541.json", "Without Skills (earlier k=3 run)"),
    ("results_2026-07-19T16-55-57-401.json", "Without Skills (k=1 run)"),
]

# (file in results/coding_agents/, condition)
CODING_FILES = [
    ("coding_agent_claude_2026-08-16T14-23-33.json", "Without Skills"),
    ("coding_agent_skills-claude_2026-08-16T11-20-54.json", "With Skills"),
    ("coding_agent_codex_2026-08-16T20-33-59.json", "Without Skills"),
    ("coding_agent_skills-codex_2026-08-16T17-09-28.json", "With Skills"),
    ("coding_agent_gemini_2026-08-17T12-23-36.json", "Without Skills"),
    ("coding_agent_skills-gemini_2026-08-17T08-54-26.json", "With Skills"),
]

DISPLAY_NAMES = {
    "gpt-5.5": "GPT-5.5",
    "gpt-5.4-mini": "GPT-5.4 Mini",
    "gemini-openrouter": "Gemini 3.5 Flash",
    "gemini-flash-lite": "Gemini 3.1 Flash-Lite",
    "claude-opus-4-8": "Claude Opus 4.8",
    "claude-sonnet-5": "Claude Sonnet 5",
    "claude-haiku-4-5": "Claude Haiku 4.5",
    "rolv": "Agentic RAG",
    "rolv-nopreexec": "Agentic RAG (no pre-exec)",
    "rolv-raw": "Claude Sonnet 4.5",
    "claude": "Claude Code",
    "codex": "Codex CLI",
    "gemini": "Gemini CLI",
}


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def fmt(value):
    return "" if value is None else f"{value:.4f}"


def attempt_success_rate(pass_k):
    """Mean share of passed attempts per prompt (older files do not store it)."""
    if pass_k.get("avgAttemptSuccessRate") is not None:
        return pass_k["avgAttemptSuccessRate"]
    rates = [sum(bool(a.get("passed")) for a in r["attempts"]) / len(r["attempts"])
             for r in pass_k["results"] if r.get("attempts")]
    return sum(rates) / len(rates) if rates else None


def grading(data, pass_k):
    """How attempts were judged in this file."""
    if data.get("grader_version"):
        return data["grader_version"]
    validated = any("validated" in a for r in pass_k["results"] for a in r.get("attempts", []))
    return "validation_v1" if validated else "execution_only"


def summarise_api():
    rows = []
    for name, condition in API_FILES:
        data = load_json(RESULTS / "api" / name)
        for model, entry in data["results"].items():
            pass_k = entry.get("passK")
            if not pass_k:
                continue
            quality = entry.get("quality", {})
            rows.append({
                "condition": condition,
                "model": model,
                "display_name": DISPLAY_NAMES.get(model, model),
                "k": data["config"]["k"],
                "n_prompts": len(pass_k["results"]),
                "pass_k": fmt(pass_k["avg"]),
                "attempt_success_rate": fmt(attempt_success_rate(pass_k)),
                "workflow": fmt(quality.get("avgWorkflow")),
                "grounding": fmt(quality.get("avgGrounding")),
                "grading": grading(data, pass_k),
                "source_file": name,
            })
    return rows


def summarise_coding_agents():
    rows = []
    for name, condition in CODING_FILES:
        data = load_json(RESULTS / "coding_agents" / name)
        attempts, metrics = data["results"], data["metrics"]
        rows.append({
            "condition": condition,
            "agent": data["agent"],
            "display_name": DISPLAY_NAMES.get(data["agent"], data["agent"]),
            "k": data["k"],
            "n_prompts": data["n_prompts"],
            "n_attempts": len(attempts),
            "pass_at_k": fmt(metrics["pass_at_k"]),
            "r_success_rate": fmt(metrics["r_success_rate"]),
            # share of attempts that passed the grader
            "completion_rate": fmt(sum(bool(a.get("grader_pass")) for a in attempts) / len(attempts)),
            "avg_grader_score": fmt(metrics["avg_score"]),
            "avg_latency_s": f"{metrics['avg_latency_ms'] / 1000:.1f}",
            "benchmark_version": data["benchmark_version"],
            "source_file": name,
        })
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", path.name, f"({len(rows)} rows)")


if __name__ == "__main__":
    write_csv(summarise_api(), RESULTS / "api_summary.csv")
    write_csv(summarise_coding_agents(), RESULTS / "coding_agents_summary.csv")

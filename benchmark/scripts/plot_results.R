#!/usr/bin/env Rscript
# plot_results.R
#
# Generates publication-ready figures from benchmark_runner.js CSV output.
# Run after benchmark_runner.js completes.
#
# Usage:
#   Rscript plot_results.R \
#     --csv results/results_<run_id>_summary.csv \
#     --out_dir results/figures/

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

option_list <- list(
  make_option("--csv",     type="character", help="Path to summary CSV from benchmark_runner"),
  make_option("--out_dir", type="character", default="results/figures")
)
opts <- parse_args(OptionParser(option_list=option_list))
dir.create(opts$out_dir, recursive=TRUE, showWarnings=FALSE)

if (is.null(opts$csv) || !file.exists(opts$csv)) {
  stop("Please provide --csv path to the summary CSV")
}

df <- read.csv(opts$csv, stringsAsFactors=FALSE)

# Clean metric names for display
metric_labels <- c(
  pass_at_k        = "Pass@k",
  workflow_success = "Workflow Success Rate",
  grounding_rate   = "Hallucination / Grounding Rate",
  routing_accuracy = "Routing Accuracy",
  routing_macro_f1 = "Routing Macro-F1",
  gene_turing_score= "Gene-Turing Score (QA)"
)
model_labels <- c(
  rolv   = "ROLV (Ours)",
  gpt4o  = "GPT-4o",
  gemini = "Gemini 1.5 Pro"
)
model_colors <- c(
  "ROLV (Ours)"     = "#2C7BB6",
  "GPT-4o"          = "#D7191C",
  "Gemini 1.5 Pro"  = "#1A9641"
)

df <- df %>%
  mutate(
    metric_label = recode(metric, !!!metric_labels),
    model_label  = recode(model,  !!!model_labels),
    value_pct    = value * 100
  )

# ── FIGURE 1: Internal Metrics Bar Chart ─────────────────────────────────────
internal_metrics <- c("pass_at_k", "workflow_success", "grounding_rate")
df_internal <- df %>% filter(metric %in% internal_metrics)

if (nrow(df_internal) > 0) {
  p1 <- ggplot(df_internal, aes(x=metric_label, y=value_pct, fill=model_label)) +
    geom_bar(stat="identity", position=position_dodge(width=0.75), width=0.65) +
    geom_text(aes(label=sprintf("%.1f%%", value_pct)),
              position=position_dodge(width=0.75), vjust=-0.4, size=3.2, fontface="bold") +
    scale_fill_manual(values=model_colors, name="Model") +
    scale_y_continuous(limits=c(0, 115), breaks=seq(0, 100, 20), labels=function(x) paste0(x, "%")) +
    labs(
      title   = "Internal Evaluation Metrics",
      subtitle = "Code quality comparison across models on Seurat PBMC workflows",
      x = NULL, y = "Score (%)"
    ) +
    theme_minimal(base_size=13) +
    theme(
      plot.title    = element_text(face="bold", size=15),
      legend.position = "top",
      axis.text.x  = element_text(size=11),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(opts$out_dir, "fig1_internal_metrics.pdf"), p1, width=9, height=6)
  ggsave(file.path(opts$out_dir, "fig1_internal_metrics.png"), p1, width=9, height=6, dpi=300)
  cat("[plot] Saved fig1_internal_metrics\n")
}

# ── FIGURE 2: Gene-Turing Score Comparison ────────────────────────────────────
df_qa <- df %>% filter(metric == "gene_turing_score")

if (nrow(df_qa) > 0) {
  p2 <- ggplot(df_qa, aes(x=reorder(model_label, value_pct), y=value_pct, fill=model_label)) +
    geom_bar(stat="identity", width=0.55) +
    geom_text(aes(label=sprintf("%.3f", value)), hjust=-0.1, size=4, fontface="bold") +
    scale_fill_manual(values=model_colors, guide="none") +
    scale_y_continuous(limits=c(0, 1.15), labels=scales::percent_format(scale=1)) +
    coord_flip() +
    labs(
      title   = "Genomics QA: Normalized Gene-Turing Score",
      subtitle = "Factual accuracy on single-cell biology questions",
      x = NULL, y = "Score"
    ) +
    theme_minimal(base_size=13) +
    theme(plot.title=element_text(face="bold", size=15))

  ggsave(file.path(opts$out_dir, "fig2_gene_turing.pdf"), p2, width=8, height=4)
  ggsave(file.path(opts$out_dir, "fig2_gene_turing.png"), p2, width=8, height=4, dpi=300)
  cat("[plot] Saved fig2_gene_turing\n")
}

# ── FIGURE 3: Routing Accuracy Metrics (ROLV only) ───────────────────────────
df_routing <- df %>% filter(metric %in% c("routing_accuracy", "routing_macro_f1") & model == "rolv")

if (nrow(df_routing) > 0) {
  p3 <- ggplot(df_routing, aes(x=metric_label, y=value_pct)) +
    geom_bar(stat="identity", fill="#2C7BB6", width=0.45) +
    geom_text(aes(label=sprintf("%.1f%%", value_pct)), vjust=-0.4, size=4.5, fontface="bold") +
    scale_y_continuous(limits=c(0,115), breaks=seq(0,100,20), labels=function(x) paste0(x,"%")) +
    labs(
      title   = "Routing Accuracy (ROLV Classifier)",
      subtitle = "Logistic regression domain classifier performance",
      x=NULL, y="Score (%)"
    ) +
    theme_minimal(base_size=13) +
    theme(
      plot.title=element_text(face="bold", size=15),
      panel.grid.major.x=element_blank()
    )

  ggsave(file.path(opts$out_dir, "fig3_routing.pdf"), p3, width=6, height=5)
  ggsave(file.path(opts$out_dir, "fig3_routing.png"), p3, width=6, height=5, dpi=300)
  cat("[plot] Saved fig3_routing\n")
}

# ── FIGURE 4: Overall Summary Radar / Heatmap ─────────────────────────────────
# Heatmap of all metrics × all models
df_all <- df %>%
  filter(metric %in% c("pass_at_k","workflow_success","grounding_rate","gene_turing_score")) %>%
  filter(model %in% c("rolv","gpt4o","gemini"))

if (nrow(df_all) > 0) {
  p4 <- ggplot(df_all, aes(x=metric_label, y=model_label, fill=value)) +
    geom_tile(color="white", linewidth=0.8) +
    geom_text(aes(label=sprintf("%.2f", value)), size=4, fontface="bold", color="white") +
    scale_fill_gradient2(low="#D7191C", mid="#FFFFBF", high="#1A9641",
                         midpoint=0.5, limits=c(0,1),
                         name="Score") +
    labs(
      title   = "Benchmark Heatmap: All Metrics × All Models",
      subtitle= "Scores normalized 0–1; green = better",
      x=NULL, y=NULL
    ) +
    theme_minimal(base_size=13) +
    theme(
      plot.title      = element_text(face="bold", size=15),
      axis.text.x     = element_text(angle=25, hjust=1),
      legend.position = "right"
    )

  ggsave(file.path(opts$out_dir, "fig4_heatmap.pdf"), p4, width=10, height=4.5)
  ggsave(file.path(opts$out_dir, "fig4_heatmap.png"), p4, width=10, height=4.5, dpi=300)
  cat("[plot] Saved fig4_heatmap\n")
}

cat("\n✅ All figures saved to:", opts$out_dir, "\n")

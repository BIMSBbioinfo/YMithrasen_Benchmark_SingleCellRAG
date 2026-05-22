setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")

tryCatch({
  suppressMessages(suppressWarnings({
# Find marker genes for all clusters with log2FC threshold of 0.25
all_markers <- FindAllMarkers(pbmc, logfc.threshold = 0.25, only.pos = TRUE)

# Filter for significant markers (adjusted p-value < 0.05)
significant_markers <- all_markers %>%
  dplyr::filter(p_val_adj < 0.05)

# Select top 5 markers per cluster based on avg_log2FC
top5_markers <- significant_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()

# Extract gene names for heatmap
top5_genes <- top5_markers$gene

# Generate heatmap of top 5 markers per cluster
DoHeatmap(pbmc, features = top5_genes) + NoLegend()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
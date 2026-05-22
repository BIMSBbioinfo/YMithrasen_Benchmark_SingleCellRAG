setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")

tryCatch({
  suppressMessages(suppressWarnings({
# Find marker genes for all clusters
pbmc.markers <- FindAllMarkers(pbmc, logfc.threshold = 0.25, min.pct = 0.25, only.pos = TRUE)

# Filter for significant markers
pbmc.markers.sig <- pbmc.markers %>%
  dplyr::filter(p_val_adj < 0.05)

# Get top 5 markers per cluster
top5_markers <- pbmc.markers.sig %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

# Extract unique gene names
genes_to_plot <- unique(top5_markers$gene)

# Scale data if not already scaled
if (length(pbmc@assays$RNA@scale.data) == 0) {
  pbmc <- ScaleData(pbmc, features = rownames(pbmc))
}

# Create heatmap
DoHeatmap(pbmc, features = genes_to_plot)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
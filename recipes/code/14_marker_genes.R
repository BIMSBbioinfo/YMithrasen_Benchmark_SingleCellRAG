library(Seurat)
library(dplyr)

# =============================================================================
# PARAMETERS
# =============================================================================
test_method     <- "wilcox"   # "wilcox" (default), "MAST", "DESeq2", "t", "LR"
logfc_threshold <- 0.25        # minimum log2FC for FindAllMarkers
min_pct         <- 0.25        # minimum fraction of cells expressing the gene
only_positive   <- TRUE        # report only upregulated markers
top_n_heatmap   <- 10          # top N markers per cluster for heatmap
output_file     <- "all_cluster_markers.csv"  # path to save marker table

# Specific cluster comparison (optional — set to NULL to skip)
cluster_ident_1 <- NULL   # e.g. 2 — cluster to test
cluster_ident_2 <- NULL   # e.g. c(0, 3) — clusters to compare against (NULL = all others)

# Genes to visualise — replace with markers relevant to your dataset
feature_genes <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                   "FCGR3A", "LYZ", "PPBP", "CD8A")
# =============================================================================

# Optional: find markers for a specific cluster comparison
if (!is.null(cluster_ident_1)) {
  cluster_markers <- FindMarkers(seu_filtered, ident.1 = cluster_ident_1,
                                 ident.2 = cluster_ident_2, test.use = test_method)
  cat("Top markers for cluster", cluster_ident_1, ":\n")
  print(head(cluster_markers, 5))
}

# Find markers for every cluster vs all others
sc_markers <- FindAllMarkers(
  seu_filtered,
  only.pos        = only_positive,
  min.pct         = min_pct,
  logfc.threshold = logfc_threshold,
  test.use        = test_method
)

# Top markers per cluster
top_markers <- sc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > logfc_threshold) %>%
  slice_head(n = 5) %>%
  ungroup()

cat("Top markers per cluster:\n")
print(top_markers)

# Save full table
write.csv(sc_markers, file = output_file, row.names = FALSE)
cat("Marker table saved to:", output_file, "\n")

# Visualisations
FeaturePlot(seu_filtered, features = feature_genes)
VlnPlot(seu_filtered,    features = feature_genes)

# Heatmap of top N markers per cluster
top_heatmap <- sc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > logfc_threshold) %>%
  slice_head(n = top_n_heatmap) %>%
  ungroup()

DoHeatmap(seu_filtered, features = top_heatmap$gene) + NoLegend()

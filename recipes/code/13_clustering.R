library(Seurat)

# =============================================================================
# PARAMETERS
# =============================================================================
resolution <- 0.5    # clustering resolution — higher = more clusters
algorithm  <- 1      # 1 = Louvain (default), 2 = Louvain with multilevel refinement,
                     # 3 = SLM, 4 = Leiden (requires leidenalg Python package)

# Suggested ranges:
#   resolution: 0.1–1.5 (start with 0.5, inspect UMAP, adjust to taste)
#   Use lower values for fewer, broader clusters
#   Use higher values for more granular clusters
# =============================================================================

seu_filtered <- FindClusters(
  seu_filtered,
  resolution = resolution,
  algorithm  = algorithm,
  verbose    = FALSE
)

cat("Cluster sizes:\n")
print(table(Idents(seu_filtered)))
cat("Total clusters:", length(levels(Idents(seu_filtered))), "\n")

DimPlot(seu_filtered, reduction = "umap", label = TRUE) +
  ggtitle(paste0("Clusters (res = ", resolution, ", algorithm = ", algorithm, ")"))

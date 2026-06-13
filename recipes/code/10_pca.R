library(Seurat)

# =============================================================================
# PARAMETERS
# =============================================================================
n_pcs_compute <- 50   # number of PCs to compute (compute more than you plan to use)
n_pcs_print   <- 5    # number of PCs to print loadings for
n_features_print <- 5 # number of top features to print per PC
# =============================================================================

seu_filtered <- RunPCA(
  seu_filtered,
  features = VariableFeatures(object = seu_filtered),
  npcs     = n_pcs_compute,
  verbose  = FALSE
)

# Print top gene loadings
print(seu_filtered[["pca"]], dims = 1:n_pcs_print, nfeatures = n_features_print)

# Visualisations
VizDimLoadings(seu_filtered, dims = 1:2, reduction = "pca")
DimPlot(seu_filtered, reduction = "pca") + NoLegend()
DimHeatmap(seu_filtered, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot — use this to choose n_pcs for downstream steps
# Pick the number of PCs where the curve flattens
ElbowPlot(seu_filtered, ndims = n_pcs_compute)
cat("Inspect the elbow plot and set n_pcs in 11_umap.R, 12_neighbours.R, and 05_doublets.R\n")

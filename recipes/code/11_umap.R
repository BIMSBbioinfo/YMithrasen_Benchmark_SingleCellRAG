library(Seurat)

# =============================================================================
# PARAMETERS
# =============================================================================
n_pcs      <- 20   # number of PCs — set based on elbow plot from 10_pca.R
reduction  <- "pca" # "pca" for standard; "harmony" after integration
# =============================================================================

seu_filtered <- RunUMAP(
  seu_filtered,
  dims      = 1:n_pcs,
  reduction = reduction,
  verbose   = FALSE
)

DimPlot(seu_filtered, reduction = "umap", label = TRUE) +
  ggtitle(paste0("UMAP (", n_pcs, " PCs from ", reduction, ")"))

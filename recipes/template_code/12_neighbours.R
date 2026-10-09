library(Seurat)

# PARAMETERS
n_pcs     <- 20      # must match value used in 11_umap.R
reduction <- "pca"   # "pca" for standard; "harmony" after Harmony integration


seu_filtered <- FindNeighbors(
  seu_filtered,
  dims      = 1:n_pcs,
  reduction = reduction,
  verbose   = FALSE
)

cat("kNN/SNN graph built using", n_pcs, "PCs from", reduction, "\n")

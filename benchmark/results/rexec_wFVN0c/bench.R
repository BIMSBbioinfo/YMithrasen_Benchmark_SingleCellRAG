setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19")

tryCatch({
  suppressMessages(suppressWarnings({
# Normalize the data
pbmc <- NormalizeData(pbmc)

# Find top 2000 variable features
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc)

# Find neighbors using PCA dimensions 1:30
pbmc <- FindNeighbors(pbmc, dims = 1:30)

# Run Leiden clustering with algorithm 4
pbmc <- FindClusters(pbmc, algorithm = 4, resolution = 0.8)

# Run UMAP using same dimensions as FindNeighbors
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualize UMAP colored by cluster
DimPlot(pbmc, group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
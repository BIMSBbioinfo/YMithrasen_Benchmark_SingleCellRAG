setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19")

tryCatch({
  suppressMessages(suppressWarnings({
# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find top 2000 variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Find neighbors using PCA
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:30)

# Run Leiden clustering
set.seed(42)
pbmc <- FindClusters(pbmc, algorithm = 4, resolution = 0.8)

# Run UMAP
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)

# Visualize UMAP colored by cluster
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
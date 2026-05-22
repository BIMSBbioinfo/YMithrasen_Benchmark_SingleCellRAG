setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19")

tryCatch({
  suppressMessages(suppressWarnings({
# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale all genes
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors using PCA
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:30)

# Run Leiden clustering
set.seed(42)
pbmc <- FindClusters(pbmc, algorithm = 4, resolution = 0.5)

# Run UMAP
set.seed(42)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
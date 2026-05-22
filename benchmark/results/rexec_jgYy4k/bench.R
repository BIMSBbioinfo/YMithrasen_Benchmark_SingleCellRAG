setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:30)

# Run Leiden clustering
pbmc <- FindClusters(pbmc, algorithm = 4, resolution = 0.8)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
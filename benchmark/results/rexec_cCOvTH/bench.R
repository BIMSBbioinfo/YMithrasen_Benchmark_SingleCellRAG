setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")

tryCatch({
  suppressMessages(suppressWarnings({
# Normalize data using LogNormalize method
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find top 2000 variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors using PCA dimensions 1:30
pbmc <- FindNeighbors(pbmc, dims = 1:30)

# Run Leiden clustering with algorithm 4 and resolution 0.8
pbmc <- FindClusters(pbmc, algorithm = 4, resolution = 0.8)

# Run UMAP using PCA dimensions 1:30
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualize UMAP colored by clusters
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
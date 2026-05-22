setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))

# SCTransform normalization on each sample independently
pbmc_sample1 <- SCTransform(pbmc_sample1, verbose = FALSE)
pbmc_sample2 <- SCTransform(pbmc_sample2, verbose = FALSE)

# Create list of samples
samples_list <- list(pbmc_sample1, pbmc_sample2)

# Select integration features
integration_features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)
samples_list <- PrepSCTIntegration(object.list = samples_list, anchor.features = integration_features)

# Find integration anchors with SCT normalization
integration_anchors <- FindIntegrationAnchors(
  object.list = samples_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  dims = 1:30
)

# Integrate data with SCT normalization
integrated_seurat <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  dims = 1:30
)

# Set default assay to integrated
DefaultAssay(integrated_seurat) <- "integrated"

# Run PCA, UMAP, and clustering on integrated data
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = FALSE)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

# Visualize UMAP split by sample
umap_split <- DimPlot(
  integrated_seurat,
  reduction = "umap",
  split.by = "orig.ident",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("UMAP Split by Sample")
print(umap_split)

# Visualize UMAP by clusters
umap_clusters <- DimPlot(
  integrated_seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.5
) + ggtitle("UMAP Clusters")
print(umap_clusters)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
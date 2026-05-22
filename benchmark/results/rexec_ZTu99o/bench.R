setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Run SCTransform on each sample independently
pbmc_sample1 <- SCTransform(pbmc_sample1, verbose = FALSE)
pbmc_sample2 <- SCTransform(pbmc_sample2, verbose = FALSE)

# Create list of samples
pbmc_list <- list(pbmc_sample1, pbmc_sample2)

# Select integration features
features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures = 3000)
pbmc_list <- PrepSCTIntegration(object.list = pbmc_list, anchor.features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = pbmc_list, normalization.method = "SCT", 
                                   anchor.features = features, verbose = FALSE)

# Integrate data
pbmc_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                                  verbose = FALSE)

# Set default assay to integrated
DefaultAssay(pbmc_integrated) <- "integrated"

# Run dimensionality reduction
pbmc_integrated <- RunPCA(pbmc_integrated, verbose = FALSE)
pbmc_integrated <- RunUMAP(pbmc_integrated, dims = 1:30, verbose = FALSE)

# Find neighbors and clusters
pbmc_integrated <- FindNeighbors(pbmc_integrated, dims = 1:30, verbose = FALSE)
pbmc_integrated <- FindClusters(pbmc_integrated, resolution = 0.5, verbose = FALSE)

# Visualize UMAP split by sample
DimPlot(pbmc_integrated, reduction = "umap", split.by = "orig.ident")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
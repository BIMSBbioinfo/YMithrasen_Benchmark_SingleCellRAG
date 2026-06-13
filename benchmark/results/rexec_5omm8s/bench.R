setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(dplyr)

# Define file paths (assuming 10X Genomics format)
file_paths <- list.files(pattern = "PBMC", full.names = TRUE)

# Read each dataset and create Seurat objects
seurat_list <- lapply(file_paths, function(path) {
  data <- Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(path))
  return(seurat_obj)
})

# Merge the three datasets
pbmc_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = c("PBMC1", "PBMC2", "PBMC3"))

# Preprocessing
pbmc_merged <- NormalizeData(pbmc_merged, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_merged <- FindVariableFeatures(pbmc_merged, selection.method = "vst", nfeatures = 2000)
pbmc_merged <- ScaleData(pbmc_merged, features = rownames(pbmc_merged))
pbmc_merged <- RunPCA(pbmc_merged, features = VariableFeatures(object = pbmc_merged))

# Determine dimensionality
ElbowPlot(pbmc_merged)

# Run UMAP and clustering
pbmc_merged <- RunUMAP(pbmc_merged, dims = 1:20)
pbmc_merged <- FindNeighbors(pbmc_merged, dims = 1:20)
pbmc_merged <- FindClusters(pbmc_merged, resolution = 0.5)

# Visualize clusters
DimPlot(pbmc_merged, reduction = "umap", label = TRUE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
```{r}
# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(SeuratData)
  library(SeuratDisk)
})

# Set seed for reproducibility
set.seed(1234)

# Define file paths
count_file <- "counts_matrix.csv"          # Replace with your actual file
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# Read raw count matrix
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC plots
pdf(file.path(output_dir, "QC_violin_plots.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data using SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Run PCA
pbmc <- RunPCA(pbmc, verbose = FALSE)

# Determine dimensionality
pdf(file.path(output_dir, "ElbowPlot.pdf"))
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)

# Save Seurat object
saveRDS(pbmc, file = file.path(output_dir, "pbmc_seurat_object.rds"))
SaveH5Seurat(pbmc, filename = file.path(output_dir, "pbmc_seurat_object.h5Seurat"))

# UMAP plot
pdf(file.path(output_dir, "UMAP_clusters.pdf"))
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of PBMCs")
dev.off()

# PCA plot
pdf(file.path(output_dir, "P
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
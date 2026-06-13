setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

# Set file paths (adjust as needed)
barcode_file <- "barcodes.tsv"
gene_file    <- "genes.tsv"
matrix_file  <- "matrix.mtx"

# Read the data
barcodes <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)$V1
genes    <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
rownames(genes) <- genes$V1
colnames(genes) <- genes$V2

# Load the sparse matrix
counts <- readMM(matrix_file)
colnames(counts) <- barcodes
rownames(counts) <- genes$V2

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization and scaling
pbmc <- SCTransform(pbmc, verbose = FALSE)

# PCA
pbmc <- RunPCA(pbmc, verbose = FALSE)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

# Save the final Seurat object
saveRDS(pbmc, file = "pbmc_final.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
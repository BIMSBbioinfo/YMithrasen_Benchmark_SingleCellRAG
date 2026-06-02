setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load necessary libraries
library(Seurat)
library(Matrix)

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create a Seurat object
pbmc <- CreateSeuratObject(counts = matrix, rownames = genes$V2, colnames = barcodes$V1)

# Perform standard quality control
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate the number of detected genes per cell and total UMI counts per cell
pbmc$detected_genes <- Matrix::rowSums(pbmc@assays$RNA@counts > 0)
pbmc$total_umi <- Matrix::rowSums(pbmc@assays$RNA@counts)

# Filter out low-quality cells
pbmc <- subset(pbmc, subset = detected_genes > 200 & detected_genes < 2500 & total_umi > 1000 & total_umi < 20000 & percent.mt < 5)

# Print the number of cells after filtering
print(paste("Number of cells after filtering: ", ncol(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
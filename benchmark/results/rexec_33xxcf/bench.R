setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required library
library(Seurat)

# Method 1: Using Read10X (recommended)
# Read10X expects a directory containing the three files
# Assuming your files are in a directory called "filtered_gene_bc_matrices" or similar
pbmc_data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc_data, 
                          project = "PBMC3k", 
                          min.cells = 3,      # Include genes detected in at least 3 cells
                          min.features = 200) # Include cells with at least 200 detected genes

# View basic information
pbmc

# Method 2: If Read10X doesn't work, manual loading
# Load the matrix manually
library(Matrix)

# Read the files
matrix_data <- readMM("path/to/matrix.mtx")
barcodes <- read.table("path/to/barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("path/to/genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix_data) <- genes$V2  # Use HGNC symbols (second column)
colnames(matrix_data) <- barcodes$V1

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix_data, 
                          project = "PBMC3k",
                          min.cells = 3,
                          min.features = 200)

# Basic QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View summary
print(paste("Number of cells:", ncol(pbmc)))
print(paste("Number of genes:", nrow(pbmc)))
head(pbmc@meta.data)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
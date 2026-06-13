setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required library
library(Seurat)

# Set the path to your data directory
# Replace with your actual path where the three files are located
data_dir <- "."

# Load the 10X data using Read10X
# This function expects barcodes.tsv, genes.tsv, and matrix.mtx in the directory
pbmc_data <- Read10X(data.dir = data_dir)

# Create a Seurat object
# min.cells: include genes detected in at least 3 cells
# min.features: include cells with at least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc_data, 
                          project = "pbmc2700", 
                          min.cells = 3, 
                          min.features = 200)

# View basic information about the object
print(pbmc)

# Optional: View the first few rows and columns of the count matrix
pbmc_data[1:10, 1:5]

# Optional: Check dimensions
cat("Number of genes:", nrow(pbmc), "\n")
cat("Number of cells:", ncol(pbmc), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
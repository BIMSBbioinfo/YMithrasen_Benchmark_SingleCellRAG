setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)

# Set working directory to where your files are located
# setwd("path/to/your/data")

# Load the sparse matrix data
matrix <- readMM("matrix.mtx")

# Load barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Load genes (feature identifiers)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# View basic information about the object
print(pbmc)

# Normalize the data using log-normalization
# scale.factor = 10000 (default)
# This performs: log1p(counts / total_counts_per_cell * 10000)
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# View normalized data (first 5 genes, first 5 cells)
pbmc[["RNA"]]@data[1:5, 1:5]

# Optional: Save the Seurat object for later use
saveRDS(pbmc, file = "pbmc_2700_normalized.rds")

# Summary statistics
cat("\nDataset Summary:\n")
cat("Number of cells:", ncol(pbmc), "\n")
cat("Number of genes:", nrow(pbmc), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
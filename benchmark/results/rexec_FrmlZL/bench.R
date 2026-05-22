setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required library
suppressMessages(suppressWarnings(library(Seurat)))

# Read 10X Genomics format data from current directory
# Read10X() automatically detects and reads barcodes.tsv, genes.tsv, and matrix.mtx
counts <- Read10X(data.dir = ".")

# Create Seurat object with quality filtering
# min.cells: Include genes detected in at least 3 cells
# min.features: Include cells with at least 200 detected genes
pbmc <- CreateSeuratObject(
  counts = counts,
  project = "PBMC2700",
  min.cells = 3,
  min.features = 200
)

# Display basic information about the Seurat object
print(pbmc)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
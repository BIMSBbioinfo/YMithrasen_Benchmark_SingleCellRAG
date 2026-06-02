setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(Matrix)

# Read in the matrix
matrix <- readMM("matrix.mtx")

# Read in the genes and barcodes
genes <- read.table("genes.tsv", colClasses = c("character", "character"), 
                     col.names = c("gene_id", "hgnc_symbol"), header = FALSE, 
                     row.names = 1, check.names = FALSE)

barcodes <- read.table("barcodes.tsv", colClasses = "character", 
                       col.names = "cell_id", header = FALSE, 
                       row.names = 1, check.names = FALSE)

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, 
                                    project = "PBMCs", 
                                    genes = genes, 
                                    cells = barcodes)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
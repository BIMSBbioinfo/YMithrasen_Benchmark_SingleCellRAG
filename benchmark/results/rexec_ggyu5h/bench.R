setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(Matrix)

# Load data
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
rownames(genes) <- genes$V2  # use HGNC symbols as row names
matrix_counts <- readMM("matrix.mtx")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix_counts, project = "PBMC", min.cells = 0, min.features = 0)
rownames(pbmc) <- genes$V2
pbmc@meta.data$barcode <- barcodes$V1

# Calculate percent mitochondrial genes
mito_genes <- grep("^MT-", rownames(pbmc), value = TRUE)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mito_genes)

# Filter cells
pbmc_filtered <- subset(pbmc, subset = nFeature_RNA >= 200 & nFeature_RNA <= 2500 & percent.mt < 5)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
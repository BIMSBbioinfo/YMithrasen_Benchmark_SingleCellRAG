setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(Matrix)
library(ggplot2)

# Load data
barcodes <- read.table("barcodes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Read sparse matrix
mat <- readMM("matrix.mtx")

# Assign row and column names
rownames(mat) <- genes$gene_symbol
colnames(mat) <- barcodes$V1

# Create Seurat object
pbmc <- CreateSeuratObject(counts = mat, project = "PBMC", min.cells = 3, min.features = 200)

# Calculate percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Feature scatter plots
p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("nCount_RNA vs percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("nCount_RNA vs nFeature_RNA")

print(p1)
print(p2)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
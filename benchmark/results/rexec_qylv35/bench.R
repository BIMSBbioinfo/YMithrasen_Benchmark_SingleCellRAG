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
matrix_file <- "matrix.mtx"
barcodes_file <- "barcodes.tsv"
genes_file <- "genes.tsv"

counts <- readMM(matrix_file)
barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- genes$V2  # use HGNC symbols as gene names
colnames(counts) <- barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Calculate percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Feature scatter plots
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  ggtitle("nCount_RNA vs percent.mt")

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  ggtitle("nCount_RNA vs nFeature_RNA")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
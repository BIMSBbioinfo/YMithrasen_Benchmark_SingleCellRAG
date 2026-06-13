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

# Calculate the number of detected genes per cell
pbmc <- CalculateFeatures(object = pbmc, features = rownames(pbmc))

# Calculate the total UMI counts per cell
pbmc <- AddModuleScore(object = pbmc, features = rownames(pbmc), name = "nUMI")

# Identify mitochondrial genes
mito.genes <- genes[grep("^MT-", genes$V2), ]$V2

# Calculate the percentage of mitochondrial genes
pbmc <- AddModuleScore(object = pbmc, features = mito.genes, name = "percent.mito")

# Filter out low-quality cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.05 & nUMI > 1000 & nUMI < 20000)

# Print the number of cells after filtering
print(paste("Number of cells after filtering: ", ncol(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
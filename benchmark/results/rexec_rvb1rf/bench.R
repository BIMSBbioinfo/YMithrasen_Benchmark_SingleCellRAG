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
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC", 
                           min.cells = 0, 
                           min.features = 0, 
                           genes = genes[, 2], 
                           cells = barcodes[, 1])

# Calculate the number of detected genes per cell
pbmc@meta.data$n_genes <- Matrix::rowSums(pbmc@assays$RNA@counts > 0)

# Calculate the total UMI counts per cell
pbmc@meta.data$n_UMI <- Matrix::rowSums(pbmc@assays$RNA@counts)

# Identify mitochondrial genes
mito.genes <- grep("^MT-", rownames(pbmc@assays$RNA@counts), value = TRUE)

# Calculate the percentage of mitochondrial genes
pbmc@meta.data$pct_mito <- (Matrix::rowSums(pbmc@assays$RNA@counts[mito.genes, ])) / 
  pbmc@meta.data$n_UMI * 100

# Filter out low-quality cells
pbmc <- subset(pbmc, subset = n_UMI > 1000 & n_UMI < 25000 & n_genes > 200 & n_genes < 2500 & pct_mito < 5)

# Plot the filtered data
plot1 <- DimPlot(pbmc, reduction = "pca")
plot2 <- VlnPlot(pbmc, features = c("n_UMI", "n_genes", "pct_mito"))
plot3 <- FeatureScatter(pbmc, feature1 = "n_UMI", feature2 = "n_genes")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
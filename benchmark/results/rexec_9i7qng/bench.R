setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================
# 1. Load the data files
# ============================================

# Read the sparse matrix in Matrix Market format
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================
# 2. Format the matrix
# ============================================

# Assign row names (genes) and column names (cells) to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================
# 3. Create Seurat object
# ============================================

# Create Seurat object with basic QC filtering
# min.cells: include genes detected in at least 3 cells
# min.features: include cells with at least 200 detected genes
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,
  min.features = 200
)

# View basic information
print(pbmc)

# ============================================
# 4. Quality control metrics (optional but recommended)
# ============================================

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Optional: Filter cells based on QC metrics
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ============================================
# 5. Normalize the data
# ============================================

# Log-normalization (default method)
# Normalizes feature expression for each cell by total expression
# Multiplies by scale factor (10,000 by default)
# Log-transforms the result
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("\nData normalized successfully!\n")

# ============================================
# 6. Identify highly variable features
# ============================================

# Find the top 2000 most variable genes
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 variable genes:", top10, "\n")

# Visualize variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================
# 7. Scale the data
# ============================================

# Scale data so that mean expression = 0 and variance = 1 across cells
# This is done for all genes by default, but you can specify features
all_genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all_genes,
  verbose = TRUE
)

cat("\nData scaled successfully!\n")
cat("Mean expression = 0, Variance = 1 across cells for all genes\n")

# ============================================
# 8. Verify scaling
# ============================================

# Check that scaling worked correctly
scaled_data <- GetAssayData(pbmc, slot = "scale.data")
cat("\nScaled data dimensions:", dim(scaled_data), "\n")

# Verify mean ≈ 0 and variance ≈ 1 for a random gene
random_gene <- sample(rownames(scaled_data), 1)
gene_values <- scaled_data[random_gene, ]
cat("Random gene:", random_gene, "\n")
cat("Mean:", round(mean(gene_values), 6), "\n")
cat("Variance:", round(var(gene_values), 6), "\n")

# ============================================
# 9. Save the processed Seurat object (optional)
# ============================================

# saveRDS(pbmc, file = "pbmc_2700_processed.rds")
# To load later: pbmc <- readRDS("pbmc_2700_processed.rds")

# ============================================
# Summary
# ============================================

cat("\n=== Processing Complete ===\n")
cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")
cat("Data normalized: LogNormalize method\n")
cat("Variable features identified:", length(VariableFeatures(pbmc)), "genes\n")
cat("Data scaled: mean = 0, variance = 1 across cells\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
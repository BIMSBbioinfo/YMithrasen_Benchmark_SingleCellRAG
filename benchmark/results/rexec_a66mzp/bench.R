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

# ============================================================================
# STEP 1: Load the data files
# ============================================================================

# Load the sparse matrix (gene expression counts)
matrix <- readMM("matrix.mtx")

# Load barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Load genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================================================
# STEP 2: Assign row and column names to the matrix
# ============================================================================

# Assign gene symbols as row names and barcodes as column names
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================================================
# STEP 3: Create Seurat object
# ============================================================================

# Create Seurat object with basic QC filtering
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("\nSeurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# STEP 4: Quality Control (Optional but recommended)
# ============================================================================

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Optional: Filter cells based on QC metrics
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ============================================================================
# STEP 5: Normalize the data
# ============================================================================

# Log-normalization (default method)
# This normalizes feature expression for each cell by total expression,
# multiplies by scale factor (10,000), and log-transforms the result
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("\nData normalized using LogNormalize method\n")

# ============================================================================
# STEP 6: Identify highly variable features (genes)
# ============================================================================

# Find the top 2000 most variable genes
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Visualize variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# STEP 7: Scale the data
# ============================================================================

# Scale data so that mean expression = 0 and variance = 1 across cells
# This is done for all genes by default, but you can specify features
all_genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all_genes,  # Scale all genes (can use VariableFeatures(pbmc) for variable genes only)
  do.center = TRUE,      # Center to mean = 0
  do.scale = TRUE        # Scale to variance = 1
)

cat("\nData scaled: mean = 0, variance = 1 across cells\n")

# ============================================================================
# STEP 8: Verify the scaling
# ============================================================================

# Check scaled data for a few genes
scaled_data <- GetAssayData(pbmc, slot = "scale.data")
cat("\nScaled data dimensions:", dim(scaled_data), "\n")

# Verify mean and variance for a sample gene
sample_gene <- rownames(scaled_data)[1]
cat("Sample gene:", sample_gene, "\n")
cat("Mean:", round(mean(scaled_data[sample_gene, ]), 6), "\n")
cat("Variance:", round(var(scaled_data[sample_gene, ]), 6), "\n")

# ============================================================================
# STEP 9: Save the processed Seurat object (optional)
# ============================================================================

# Save for future use
saveRDS(pbmc, file = "pbmc_processed.rds")
cat("\nProcessed Seurat object saved as 'pbmc_processed.rds'\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Processing Complete ===\n")
cat("Final object contains:\n")
cat("- Cells:", ncol(pbmc), "\n")
cat("- Genes:", nrow(pbmc), "\n")
cat("- Variable features:", length(VariableFeatures(pbmc)), "\n")
cat("- Scaled features:", nrow(GetAssayData(pbmc, slot = "scale.data")), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
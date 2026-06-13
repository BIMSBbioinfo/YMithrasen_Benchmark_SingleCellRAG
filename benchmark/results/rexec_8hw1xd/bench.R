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
library(ggplot2)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# 1. Load the data files
# ============================================================================

# Read the sparse matrix (Matrix Market format)
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================================================
# 2. Assign row and column names to the matrix
# ============================================================================

# Assign gene symbols as row names and barcodes as column names
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix)[1], "genes x", dim(matrix)[2], "cells\n")

# ============================================================================
# 3. Create Seurat object
# ============================================================================

# Create Seurat object with basic QC filtering
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# Display basic information
print(pbmc)

# ============================================================================
# 4. Basic QC metrics (optional but recommended)
# ============================================================================

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ============================================================================
# 5. Normalize the data
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 6. Identify highly variable features using VST method
# ============================================================================

pbmc <- FindVariableFeatures(
  pbmc, 
  selection.method = "vst",  # Variance Stabilizing Transformation
  nfeatures = 2000           # Select top 2000 features
)

# Get the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 highly variable genes:\n")
print(top10)

# ============================================================================
# 7. Visualize variable features
# ============================================================================

# Create the variable feature plot with top 10 labeled
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(
  plot = plot1, 
  points = top10, 
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

# Display the plot
print(plot2)

# Save the plot (optional)
ggsave("variable_features_plot.png", plot = plot2, width = 10, height = 6, dpi = 300)

# ============================================================================
# 8. Additional visualization: Heatmap of top variable genes (optional)
# ============================================================================

# Get top 20 variable genes for heatmap
top20 <- head(VariableFeatures(pbmc), 20)

# Create a simple heatmap showing expression of top variable genes
DoHeatmap(
  subset(pbmc, downsample = 100),  # Downsample to 100 cells for clarity
  features = top20,
  size = 3
) + NoLegend()

# ============================================================================
# 9. Summary statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Total cells after filtering:", ncol(pbmc), "\n")
cat("Total genes after filtering:", nrow(pbmc), "\n")
cat("Number of highly variable features:", length(VariableFeatures(pbmc)), "\n")

# Get variable feature information
variable_features_info <- HVFInfo(pbmc)
cat("\nMean variance of top 2000 features:", 
    mean(variable_features_info[VariableFeatures(pbmc), "variance.standardized"]), "\n")

# ============================================================================
# 10. Save the Seurat object (optional)
# ============================================================================

saveRDS(pbmc, file = "pbmc_2700_seurat.rds")
cat("\nSeurat object saved as 'pbmc_2700_seurat.rds'\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
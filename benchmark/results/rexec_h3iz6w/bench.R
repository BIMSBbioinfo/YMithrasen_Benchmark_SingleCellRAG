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

# ============================================================
# STEP 1: Load the data files
# ============================================================

# Load the sparse matrix (UMI counts)
matrix <- readMM("matrix.mtx")

# Load barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Load genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================================
# STEP 2: Create Seurat object
# ============================================================

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("\nSeurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================
# STEP 3: Quality Control and Filtering
# ============================================================

# Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================
# STEP 4: Normalization
# ============================================================

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("Data normalized\n")

# ============================================================
# STEP 5: Identify highly variable features (genes)
# ============================================================

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================
# STEP 6: Scaling the data
# ============================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaled\n")

# ============================================================
# STEP 7: Perform PCA on top 2000 variable genes
# ============================================================

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = 50,  # Compute 50 principal components
  verbose = FALSE
)

cat("\nPCA completed\n")

# Print the top genes associated with the first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA results
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# PCA plot colored by total UMI counts
DimPlot(pbmc, reduction = "pca")

# ============================================================
# STEP 8: Elbow Plot to determine number of PCs
# ============================================================

# Create elbow plot
elbow_plot <- ElbowPlot(pbmc, ndims = 50)
print(elbow_plot)

# Save the elbow plot
ggsave("elbow_plot.png", plot = elbow_plot, width = 8, height = 6, dpi = 300)

# ============================================================
# STEP 9: Additional PCA visualizations
# ============================================================

# Heatmap of top PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ============================================================
# STEP 10: Determine dimensionality (optional automated approach)
# ============================================================

# Calculate standard deviation for each PC
pca_stdev <- pbmc@reductions$pca@stdev

# Calculate percentage of variance explained
pct <- pca_stdev / sum(pca_stdev) * 100

# Calculate cumulative variance
cumulative <- cumsum(pct)

# Find PC where cumulative variance reaches 90%
pc_90 <- which(cumulative > 90)[1]

# Find PC where change in variance is less than 0.1%
co1 <- which(diff(pct) < 0.1)[1]

# Find PC where cumulative change is less than 0.1%
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1

# Minimum of the three methods
pcs_to_use <- min(pc_90, co1, co2)

cat("\nSuggested number of PCs to use:", pcs_to_use, "\n")
cat("PC at 90% variance:", pc_90, "\n")
cat("PC where variance change < 0.1%:", co1, "\n")

# ============================================================
# STEP 11: Save the Seurat object
# ============================================================

saveRDS(pbmc, file = "pbmc_2700_processed.rds")
cat("\nSeurat object saved as 'pbmc_2700_processed.rds'\n")

# ============================================================
# Summary Statistics
# ============================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Variable features identified:", length(VariableFeatures(pbmc)), "\n")
cat("PCs computed:", 50, "\n")
cat("Recommended PCs for downstream analysis:", pcs_to_use, "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
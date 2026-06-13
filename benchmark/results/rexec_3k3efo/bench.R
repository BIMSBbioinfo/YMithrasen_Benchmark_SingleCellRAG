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
library(dplyr)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data files
# ============================================================================

# Load the sparse matrix (gene expression data)
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

# ============================================================================
# STEP 2: Create Seurat object
# ============================================================================

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("\nSeurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# STEP 3: Quality Control and Filtering
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# View QC statistics
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

cat("\nAfter filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 4: Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

cat("Data normalized\n")

# ============================================================================
# STEP 5: Identify highly variable features (genes)
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# STEP 6: Scaling the data
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaled\n")

# ============================================================================
# STEP 7: Perform PCA on the top 2000 variable genes
# ============================================================================

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50,  # Compute 50 principal components
               verbose = FALSE)

cat("\nPCA completed\n")

# Print the top genes associated with the first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA results
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# PCA plot colored by gene expression
DimPlot(pbmc, reduction = "pca")

# Heatmap of top PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ============================================================================
# STEP 8: Elbow Plot - Determine number of PCs to use
# ============================================================================

# Create elbow plot
elbow_plot <- ElbowPlot(pbmc, ndims = 50)
print(elbow_plot)

# Enhanced elbow plot with additional information
p <- ElbowPlot(pbmc, ndims = 50) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "blue") +
  labs(title = "PCA Elbow Plot for PBMC Data",
       subtitle = "Red line: PC=10, Blue line: PC=15",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme_minimal()

print(p)

# Calculate percentage of variance explained
pca_variance <- pbmc@reductions$pca@stdev^2
variance_explained <- pca_variance / sum(pca_variance) * 100
cumulative_variance <- cumsum(variance_explained)

# Print variance statistics
cat("\nVariance explained by first 20 PCs:\n")
variance_df <- data.frame(
  PC = 1:20,
  Variance = round(variance_explained[1:20], 2),
  Cumulative = round(cumulative_variance[1:20], 2)
)
print(variance_df)

# Plot cumulative variance
cumvar_plot <- ggplot(data.frame(PC = 1:50, 
                                  CumVar = cumulative_variance[1:50]), 
                      aes(x = PC, y = CumVar)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
  labs(title = "Cumulative Variance Explained",
       x = "Principal Component",
       y = "Cumulative Variance (%)") +
  theme_minimal()

print(cumvar_plot)

# ============================================================================
# STEP 9: Summary and Recommendations
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Variable features used:", length(VariableFeatures(pbmc)), "\n")
cat("PCs computed:", 50, "\n")
cat("\nVariance explained by PC10:", round(cumulative_variance[10], 2), "%\n")
cat("Variance explained by PC15:", round(cumulative_variance[15], 2), "%\n")
cat("Variance explained by PC20:", round(cumulative_variance[20], 2), "%\n")

cat("\nRECOMMENDATION: Based on the elbow plot, consider using 10-15 PCs for downstream analysis\n")

# Save the Seurat object for future use
saveRDS(pbmc, file = "pbmc_2700_processed.rds")
cat("\nSeurat object saved as 'pbmc_2700_processed.rds'\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
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
# STEP 9: Additional visualizations for PC selection
# ============================================================

# Heatmap of top genes for first 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# JackStraw plot (optional - computationally intensive)
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)

# ============================================================
# STEP 10: Summary statistics
# ============================================================

cat("\n=== Analysis Summary ===\n")
cat("Final number of cells:", ncol(pbmc), "\n")
cat("Final number of genes:", nrow(pbmc), "\n")
cat("Number of variable features:", length(VariableFeatures(pbmc)), "\n")
cat("Number of PCs computed:", ncol(Embeddings(pbmc, "pca")), "\n")

# Calculate variance explained by each PC
pca_variance <- pbmc@reductions$pca@stdev^2
pca_variance_pct <- (pca_variance / sum(pca_variance)) * 100
cumulative_variance <- cumsum(pca_variance_pct)

cat("\nVariance explained by first 10 PCs:\n")
for(i in 1:10) {
  cat(sprintf("PC%d: %.2f%% (Cumulative: %.2f%%)\n", 
              i, pca_variance_pct[i], cumulative_variance[i]))
}

# Save the Seurat object for future use
saveRDS(pbmc, file = "pbmc2700_processed.rds")
cat("\nSeurat object saved as 'pbmc2700_processed.rds'\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
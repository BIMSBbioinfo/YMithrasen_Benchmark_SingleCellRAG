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

# Create Seurat object with basic filtering
# min.cells: include genes detected in at least 3 cells
# min.features: include cells with at least 200 detected genes
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,
  min.features = 200
)

cat("\nSeurat object created successfully!\n")
print(pbmc)

# ============================================================================
# STEP 3: Quality Control (QC) Metrics
# ============================================================================

# Calculate mitochondrial percentage
# Mitochondrial genes start with "MT-"
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1)

# Scatter plots to visualize feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# ============================================================================
# STEP 4: Filter cells based on QC metrics
# ============================================================================

# Filter cells:
# - Keep cells with feature count between 200 and 2500
# - Keep cells with mitochondrial percentage < 5%
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)

cat("\nAfter filtering:\n")
print(pbmc)

# ============================================================================
# STEP 5: Normalization
# ============================================================================

# Normalize the data using LogNormalize method
# This normalizes feature expression for each cell by total expression,
# multiplies by scale factor (10,000), and log-transforms
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

cat("\nData normalized successfully!\n")

# ============================================================================
# STEP 6: Identify highly variable features (genes)
# ============================================================================

# Find the top 2000 most variable genes
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 variable genes:\n")
print(top10)

# Plot variable features with labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, 
                     points = top10, 
                     repel = TRUE,
                     xnudge = 0,
                     ynudge = 0)
print(plot2)

# ============================================================================
# STEP 7: Scale the data
# ============================================================================

# Scale data (z-score transformation)
# This step is required before PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("\nData scaled successfully!\n")

# ============================================================================
# STEP 8: Perform PCA
# ============================================================================

# Run PCA on the top 2000 variable genes
pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50,  # Compute 50 principal components
               verbose = TRUE)

# Print the top genes associated with the first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# ============================================================================
# STEP 9: Visualize PCA results
# ============================================================================

# Visualize loadings of top genes for first two PCs
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# PCA plot (cells colored by PC scores)
DimPlot(pbmc, reduction = "pca")

# Heatmap of top genes for first 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ============================================================================
# STEP 10: Elbow Plot - Determine number of PCs
# ============================================================================

# Create elbow plot to visualize standard deviation explained by each PC
elbow_plot <- ElbowPlot(pbmc, ndims = 50)
print(elbow_plot)

# Enhanced elbow plot with more details
elbow_data <- pbmc@reductions$pca@stdev
elbow_df <- data.frame(
  PC = 1:length(elbow_data),
  StdDev = elbow_data,
  VarianceExplained = (elbow_data^2) / sum(elbow_data^2) * 100
)

# Plot with variance explained
p1 <- ggplot(elbow_df[1:30,], aes(x = PC, y = StdDev)) +
  geom_point(size = 2, color = "steelblue") +
  geom_line(color = "steelblue") +
  theme_minimal() +
  labs(title = "Elbow Plot - Standard Deviation",
       x = "Principal Component",
       y = "Standard Deviation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(elbow_df[1:30,], aes(x = PC, y = cumsum(VarianceExplained))) +
  geom_point(size = 2, color = "darkred") +
  geom_line(color = "darkred") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "Cumulative Variance Explained",
       x = "Principal Component",
       y = "Cumulative Variance (%)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)
print(p2)

# ============================================================================
# STEP 11: Summary statistics
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells after QC:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Variable features used:", length(VariableFeatures(pbmc)), "\n")
cat("PCs computed:", ncol(pbmc@reductions$pca), "\n")
cat("\nVariance explained by first 10 PCs:\n")
print(round(elbow_df$VarianceExplained[1:10], 2))
cat("\nCumulative variance by PC 10:", 
    round(sum(elbow_df$VarianceExplained[1:10]), 2), "%\n")
cat("Cumulative variance by PC 20:", 
    round(sum(elbow_df$VarianceExplained[1:20]), 2), "%\n")

# ============================================================================
# STEP 12: Save the Seurat object (optional)
# ============================================================================

# Save the processed Seurat object for future use
# saveRDS(pbmc, file = "pbmc_2700_processed.rds")
# To load later: pbmc <- readRDS("pbmc_2700_processed.rds")

cat("\n=== Analysis Complete! ===\n")
cat("Based on the elbow plot, consider using 10-20 PCs for downstream analysis.\n")
cat("Look for the 'elbow' where the curve flattens.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
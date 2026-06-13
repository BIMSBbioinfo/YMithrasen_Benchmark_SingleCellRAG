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
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Display data dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================================================
# 2. CREATE SEURAT OBJECT
# ============================================================================

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("\nSeurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 3. QUALITY CONTROL (Optional but recommended)
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells (adjust thresholds based on your data)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 4. NORMALIZATION AND SCALING
# ============================================================================

# Normalize the data (log-normalization)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("\nData normalized and scaled\n")

# ============================================================================
# 5. DIMENSIONALITY REDUCTION - PCA
# ============================================================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Elbow plot to determine number of PCs
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 6. CLUSTERING
# ============================================================================

# Build KNN graph using 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Cluster cells with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("\nClustering complete. Number of clusters:", length(unique(Idents(pbmc))), "\n")

# View cluster assignments
table(Idents(pbmc))

# ============================================================================
# 7. UMAP VISUALIZATION
# ============================================================================

# Run UMAP using 10 PCs
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Create UMAP plot
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6, pt.size = 1) +
  ggtitle("PBMC Clusters - UMAP") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right"
  )

print(p1)

# Save the plot
ggsave("pbmc_umap_clusters.png", plot = p1, width = 10, height = 8, dpi = 300)

# ============================================================================
# 8. ADDITIONAL VISUALIZATIONS
# ============================================================================

# UMAP colored by number of genes detected
p2 <- FeaturePlot(pbmc, features = "nFeature_RNA", reduction = "umap") +
  ggtitle("Number of Genes Detected")
print(p2)

# UMAP colored by UMI counts
p3 <- FeaturePlot(pbmc, features = "nCount_RNA", reduction = "umap") +
  ggtitle("UMI Counts")
print(p3)

# ============================================================================
# 9. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analyzed.rds")

# Export cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")

# ============================================================================
# 10. SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("PCs used for clustering:", 10, "\n")
cat("Clustering resolution:", 0.5, "\n")
cat("\nCells per cluster:\n")
print(table(Idents(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
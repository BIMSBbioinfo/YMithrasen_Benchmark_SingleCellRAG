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
library(dplyr)
library(ggplot2)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (feature identifiers)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_name")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_name
colnames(matrix) <- barcodes$barcode

# Display data dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================================================
# 2. CREATE SEURAT OBJECT
# ============================================================================

# Create Seurat object with basic filtering
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 3. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature scatter plots for QC
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 4. NORMALIZATION
# ============================================================================

# Normalize the data using log-normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 5. FEATURE SELECTION
# ============================================================================

# Identify highly variable features (genes)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 6. SCALING
# ============================================================================

# Scale the data (linear transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 7. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 8. CLUSTERING
# ============================================================================

# Find neighbors using the first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check cluster assignments
cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")
table(Idents(pbmc))

# ============================================================================
# 9. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP)
# ============================================================================

# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters on UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc2700_clustered.rds")

# Export cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

# ============================================================================
# 11. OPTIONAL: FIND CLUSTER MARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# Visualize top marker genes
VlnPlot(pbmc, features = c("MS4A1", "CD79A", "CD8A", "CD3D"))
FeaturePlot(pbmc, features = c("MS4A1", "CD79A", "CD8A", "CD3D"))

# Heatmap of top markers
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()

cat("\nAnalysis complete!\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
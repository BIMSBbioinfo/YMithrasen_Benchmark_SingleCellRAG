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
library(dplyr)
library(ggplot2)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data
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

# ============================================================================
# STEP 2: Create Seurat object
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
# STEP 3: Quality Control (QC)
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
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 4: Normalization
# ============================================================================

# Normalize the data using log-normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# STEP 5: Feature Selection (Identify Highly Variable Features)
# ============================================================================

# Identify the 2000 most variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# STEP 6: Scaling
# ============================================================================

# Scale the data (z-score transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 7: Linear Dimensional Reduction (PCA)
# ============================================================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine dimensionality (Elbow plot)
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# STEP 8: Clustering
# ============================================================================

# Construct k-nearest neighbor graph
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters using Louvain algorithm with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check cluster assignments
cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")
table(Idents(pbmc))

# ============================================================================
# STEP 9: Non-linear Dimensional Reduction (UMAP/tSNE)
# ============================================================================

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters on UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Clusters (Resolution 0.5)")

# Run tSNE (optional alternative visualization)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5)

# ============================================================================
# STEP 10: Find Cluster Biomarkers (Optional but recommended)
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# Visualize top marker genes
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5$gene) + NoLegend()

# ============================================================================
# STEP 11: Visualize Known Cell Type Markers
# ============================================================================

# Common PBMC markers
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Violin plots for specific markers
VlnPlot(pbmc, features = c("MS4A1", "CD3E", "CD14"))

# Dot plot for marker expression across clusters
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")) + 
  RotatedAxis()

# ============================================================================
# STEP 12: Save Results
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc2700_analyzed.rds")

# Export cluster assignments
write.csv(data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
), file = "cluster_assignments.csv", row.names = FALSE)

# Export marker genes
write.csv(pbmc.markers, file = "cluster_markers.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
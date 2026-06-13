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

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes (cell IDs)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene IDs and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

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

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("QC_violin_plots.pdf", width = 12, height = 4)

# Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("QC_scatter_plots.pdf", width = 10, height = 4)

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

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
ggsave("variable_features.pdf", width = 10, height = 6)

# ============================================================================
# 6. SCALING
# ============================================================================

# Scale the data (z-score transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 7. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# PCA plots
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("PCA_loadings.pdf", width = 10, height = 8)

DimPlot(pbmc, reduction = "pca")
ggsave("PCA_plot.pdf", width = 8, height = 6)

# Heatmap of PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("PCA_heatmap.pdf", width = 12, height = 15)

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)
ggsave("elbow_plot.pdf", width = 8, height = 6)

# ============================================================================
# 8. CLUSTERING
# ============================================================================

# Find neighbors using first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check cluster assignments
cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")
table(Idents(pbmc))

# ============================================================================
# 9. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP)
# ============================================================================

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP with clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("UMAP_clusters.pdf", width = 10, height = 8)

# ============================================================================
# 10. FIND CLUSTER MARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# Save marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "top5_markers_per_cluster.csv", row.names = FALSE)

# ============================================================================
# 11. VISUALIZATION OF KEY MARKERS
# ============================================================================

# Visualize canonical PBMC markers
canonical_markers <- c("CD3D", "CD8A", "CD4",      # T cells
                       "CD14", "LYZ",               # Monocytes
                       "MS4A1", "CD79A",            # B cells
                       "GNLY", "NKG7",              # NK cells
                       "FCGR3A", "MS4A7",           # FCGR3A+ Monocytes
                       "FCER1A", "CST3")            # Dendritic cells

# Feature plots
FeaturePlot(pbmc, features = canonical_markers[1:9], ncol = 3)
ggsave("feature_plot_markers.pdf", width = 15, height = 15)

# Dot plot
DotPlot(pbmc, features = canonical_markers) + RotatedAxis()
ggsave("dotplot_markers.pdf", width = 12, height = 6)

# Violin plot for top markers
VlnPlot(pbmc, features = canonical_markers[1:6], ncol = 3)
ggsave("violin_plot_markers.pdf", width = 15, height = 8)

# Heatmap of top markers
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()
ggsave("heatmap_top10_markers.pdf", width = 12, height = 15)

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# Save cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

# Summary statistics
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(pbmc)))

cat("\nAnalysis complete! All results saved.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
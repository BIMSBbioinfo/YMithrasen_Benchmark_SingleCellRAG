setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# ============================================================================
# 1. DATA LOADING
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes and genes
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes$V2  # Use HGNC gene symbols
colnames(matrix) <- barcodes$V1  # Use cell barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

print(paste("Loaded", ncol(pbmc), "cells and", nrow(pbmc), "genes"))

# ============================================================================
# 2. QUALITY CONTROL
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

print(paste("After filtering:", ncol(pbmc), "cells remaining"))

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
print("Top 10 variable genes:")
print(top10)

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
ggsave("variable_features.pdf", width = 10, height = 6)

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("PCA_loadings.pdf", width = 10, height = 8)

DimPlot(pbmc, reduction = "pca")
ggsave("PCA_plot.pdf", width = 8, height = 6)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ggsave("PCA_heatmap.pdf", width = 12, height = 15)

# Determine dimensionality with elbow plot
ElbowPlot(pbmc, ndims = 50)
ggsave("elbow_plot.pdf", width = 8, height = 6)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs
head(Idents(pbmc), 5)
table(Idents(pbmc))

# ============================================================================
# 8. UMAP VISUALIZATION
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("UMAP_clusters.pdf", width = 10, height = 8)

# ============================================================================
# 9. FIND ALL MARKER GENES
# ============================================================================

# Find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# View top markers
print("Top markers per cluster:")
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# Save all markers to file
write.csv(pbmc.markers, "all_cluster_markers.csv", row.names = FALSE)

# Get top 10 markers per cluster
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

write.csv(top10_markers, "top10_markers_per_cluster.csv", row.names = FALSE)

# ============================================================================
# 10. VISUALIZE MARKER GENES
# ============================================================================

# Get top 2 markers per cluster for visualization
top2_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)

# Feature plots for top markers
FeaturePlot(pbmc, features = head(top2_markers$gene, 9))
ggsave("top_markers_featureplot.pdf", width = 15, height = 15)

# Violin plots for top markers
VlnPlot(pbmc, features = head(top2_markers$gene, 9))
ggsave("top_markers_violinplot.pdf", width = 15, height = 15)

# Heatmap of top 10 markers per cluster
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()
ggsave("markers_heatmap.pdf", width = 12, height = 15)

# Dot plot of top 5 markers per cluster
top5_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)

DotPlot(pbmc, features = unique(top5_markers$gene)) + 
    RotatedAxis()
ggsave("markers_dotplot.pdf", width = 15, height = 6)

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# Print summary statistics
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(pbmc)))
cat("\nTotal marker genes found:", nrow(pbmc.markers), "\n")
cat("Marker genes per cluster:\n")
print(table(pbmc.markers$cluster))

cat("\nAnalysis complete! Check the generated PDF files for visualizations.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
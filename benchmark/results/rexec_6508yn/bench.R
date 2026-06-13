setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# ============================================================================
# PBMC Single-Cell RNA-seq Analysis Pipeline
# ============================================================================
# This script performs a complete analysis of 2,700 PBMCs including:
# - Data loading from Matrix Market format
# - Quality control and filtering
# - Normalization and scaling
# - Dimensionality reduction (PCA)
# - Clustering
# - UMAP visualization
# - Saving final Seurat object
# ============================================================================

# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading data from Matrix Market format...\n")

# Read the sparse matrix (genes x cells)
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

cat(sprintf("Loaded matrix: %d genes x %d cells\n", nrow(matrix), ncol(matrix)))

# ============================================================================
# 2. CREATE SEURAT OBJECT
# ============================================================================

cat("\nCreating Seurat object...\n")

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Keep genes expressed in at least 3 cells
  min.features = 200  # Keep cells with at least 200 detected genes
)

cat(sprintf("Seurat object created: %d genes x %d cells\n", 
            nrow(pbmc), ncol(pbmc)))

# ============================================================================
# 3. QUALITY CONTROL
# ============================================================================

cat("\nPerforming quality control...\n")

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics before filtering
cat("\nQC Metrics Summary (before filtering):\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Create QC plots
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
ggsave("qc_violin_before_filtering.pdf", p1, width = 12, height = 4)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 2500, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 200, linetype = "dashed", color = "red")
ggsave("qc_scatter_before_filtering.pdf", p2, width = 6, height = 5)

# ============================================================================
# 4. FILTER LOW-QUALITY CELLS
# ============================================================================

cat("\nFiltering low-quality cells...\n")

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = 
                 nFeature_RNA > 200 &      # Minimum genes per cell
                 nFeature_RNA < 2500 &     # Maximum genes per cell (doublet filter)
                 percent.mt < 5)           # Maximum mitochondrial percentage

cat(sprintf("After filtering: %d genes x %d cells\n", 
            nrow(pbmc), ncol(pbmc)))

# ============================================================================
# 5. NORMALIZATION
# ============================================================================

cat("\nNormalizing data...\n")

# Log-normalize the data
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

# ============================================================================
# 6. FEATURE SELECTION
# ============================================================================

cat("\nIdentifying highly variable features...\n")

# Find highly variable genes
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst",
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
p3 <- VariableFeaturePlot(pbmc)
p3 <- LabelPoints(plot = p3, points = top10, repel = TRUE)
ggsave("variable_features.pdf", p3, width = 8, height = 6)

# ============================================================================
# 7. SCALING
# ============================================================================

cat("\nScaling data...\n")

# Scale the data (z-score transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 8. DIMENSIONALITY REDUCTION - PCA
# ============================================================================

cat("\nPerforming PCA...\n")

# Run PCA
pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50,
               verbose = FALSE)

# Visualize PCA results
p4 <- DimPlot(pbmc, reduction = "pca")
ggsave("pca_plot.pdf", p4, width = 7, height = 6)

# Elbow plot to determine dimensionality
p5 <- ElbowPlot(pbmc, ndims = 50)
ggsave("elbow_plot.pdf", p5, width = 7, height = 5)

# Print top genes for first 5 PCs
cat("\nTop genes for first 5 PCs:\n")
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# ============================================================================
# 9. DETERMINE DIMENSIONALITY
# ============================================================================

# Use first 30 PCs for downstream analysis (adjust based on elbow plot)
dims_to_use <- 1:30
cat(sprintf("\nUsing %d principal components for clustering\n", length(dims_to_use)))

# ============================================================================
# 10. CLUSTERING
# ============================================================================

cat("\nPerforming clustering...\n")

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = dims_to_use)

# Find clusters (Louvain algorithm)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Print cluster sizes
cat("\nCluster sizes:\n")
print(table(Idents(pbmc)))

# ============================================================================
# 11. UMAP VISUALIZATION
# ============================================================================

cat("\nComputing UMAP...\n")

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = dims_to_use, verbose = FALSE)

# Visualize UMAP
p6 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("PBMC Clusters (UMAP)")
ggsave("umap_clusters.pdf", p6, width = 8, height = 7)

# UMAP colored by QC metrics
p7 <- FeaturePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  reduction = "umap", ncol = 3)
ggsave("umap_qc_metrics.pdf", p7, width = 15, height = 5)

# ============================================================================
# 12. FIND CLUSTER MARKERS (Optional but recommended)
# ============================================================================

cat("\nFinding cluster markers...\n")

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25,
                               verbose = FALSE)

# Get top 5 markers per cluster
top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

cat("\nTop 5 markers per cluster:\n")
print(top5_markers)

# Save marker genes to file
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)
write.csv(top5_markers, "top5_markers_per_cluster.csv", row.names = FALSE)

# Visualize top markers
if (nrow(top5_markers) > 0) {
  p8 <- DoHeatmap(pbmc, features = top5_markers$gene, size = 3) +
    NoLegend()
  ggsave("marker_heatmap.pdf", p8, width = 12, height = 10)
}

# ============================================================================
# 13. SAVE FINAL SEURAT OBJECT
# ============================================================================

cat("\nSaving final Seurat object...\n")

# Save the complete Seurat object
saveRDS(pbmc, file = "pbmc_final.rds")

cat("\n✓ Analysis complete! Seurat object saved as 'pbmc_final.rds'\n")

# ============================================================================
# 14. SUMMARY REPORT
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Final dataset: %d genes x %d cells\n", nrow(pbmc), ncol(pbmc)))
cat(sprintf("Number of clusters: %d\n", length(unique(Idents(pbmc)))))
cat(sprintf("PCs used: %d\n", length(dims_to_use)))
cat(sprintf("Resolution: 0.5\n"))
cat("\nOutput files generated:\n")
cat("  - pbmc_final.rds (Seurat object)\n")
cat("  - qc_violin_before_filtering.pdf\n")
cat("  - qc_scatter_before_filtering.pdf\n")
cat("  - variable_features.pdf\n")
cat("  - pca_plot.pdf\n")
cat("  - elbow_plot.pdf\n")
cat("  - umap_clusters.pdf\n")
cat("  - umap_qc_metrics.pdf\n")
cat("  - marker_heatmap.pdf\n")
cat("  - cluster_markers.csv\n")
cat("  - top5_markers_per_cluster.csv\n")
cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# OPTIONAL: Load the saved object later
# ============================================================================
# To load the saved object in a new session:
# pbmc <- readRDS("pbmc_final.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
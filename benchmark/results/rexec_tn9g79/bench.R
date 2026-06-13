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
# STEP 1: Load Data
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
# STEP 2: Create Seurat Object
# ============================================================================
cat("\nCreating Seurat object...\n")

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Keep genes expressed in at least 3 cells
  min.features = 200  # Keep cells with at least 200 detected genes
)

cat(sprintf("Seurat object created: %d genes x %d cells\n", 
            nrow(pbmc), ncol(pbmc)))

# ============================================================================
# STEP 3: Quality Control Metrics
# ============================================================================
cat("\nCalculating QC metrics...\n")

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate percentage of ribosomal genes
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# Visualize QC metrics before filtering
cat("\nQC Metrics Summary (before filtering):\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Create QC plots
pdf("qc_plots_before_filtering.pdf", width = 12, height = 4)
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,
        pt.size = 0.1)
dev.off()

# Feature-feature relationships
pdf("qc_scatter_plots.pdf", width = 10, height = 4)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# ============================================================================
# STEP 4: Filter Low-Quality Cells
# ============================================================================
cat("\nFiltering low-quality cells...\n")

# Filter cells based on QC metrics
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)

cat(sprintf("After filtering: %d cells retained\n", ncol(pbmc)))

# ============================================================================
# STEP 5: Normalization
# ============================================================================
cat("\nNormalizing data...\n")

# Log-normalization (default method)
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# ============================================================================
# STEP 6: Feature Selection
# ============================================================================
cat("\nIdentifying highly variable features...\n")

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
pdf("variable_features.pdf", width = 10, height = 5)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()

# ============================================================================
# STEP 7: Scaling
# ============================================================================
cat("\nScaling data...\n")

# Scale all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 8: Dimensionality Reduction - PCA
# ============================================================================
cat("\nPerforming PCA...\n")

pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50,
               verbose = FALSE)

# Visualize PCA results
pdf("pca_plots.pdf", width = 12, height = 8)
print(DimPlot(pbmc, reduction = "pca"))
print(DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE))
dev.off()

# Elbow plot to determine dimensionality
pdf("elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(pbmc, ndims = 50)
dev.off()

# ============================================================================
# STEP 9: Clustering
# ============================================================================
cat("\nPerforming clustering...\n")

# Determine the number of PCs to use (typically 10-30 for PBMCs)
n_pcs <- 20

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)

# Find clusters (resolution parameter affects granularity)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat(sprintf("Identified %d clusters\n", 
            length(unique(Idents(pbmc)))))

# ============================================================================
# STEP 10: UMAP Visualization
# ============================================================================
cat("\nComputing UMAP...\n")

pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

# Create UMAP plots
pdf("umap_clusters.pdf", width = 10, height = 8)
DimPlot(pbmc, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5) + 
  ggtitle("PBMC Clusters - UMAP")
dev.off()

# UMAP colored by QC metrics
pdf("umap_qc_metrics.pdf", width = 15, height = 5)
p1 <- FeaturePlot(pbmc, features = "nFeature_RNA", reduction = "umap")
p2 <- FeaturePlot(pbmc, features = "nCount_RNA", reduction = "umap")
p3 <- FeaturePlot(pbmc, features = "percent.mt", reduction = "umap")
print(p1 + p2 + p3)
dev.off()

# ============================================================================
# STEP 11: Find Cluster Markers (Optional but Recommended)
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

# Save markers to CSV
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)
write.csv(top5_markers, "top5_markers_per_cluster.csv", row.names = FALSE)

# Visualize top markers
pdf("top_markers_heatmap.pdf", width = 12, height = 10)
DoHeatmap(pbmc, features = top5_markers$gene) + 
  theme(text = element_text(size = 8))
dev.off()

# ============================================================================
# STEP 12: Visualize Known Cell Type Markers
# ============================================================================
cat("\nVisualizing known cell type markers...\n")

# Common PBMC markers
markers_to_plot <- c(
  "CD3D", "CD3E",      # T cells
  "CD4",               # CD4+ T cells
  "CD8A",              # CD8+ T cells
  "CD79A", "MS4A1",    # B cells
  "CD14", "LYZ",       # Monocytes
  "FCGR3A", "MS4A7",   # FCGR3A+ Monocytes
  "GNLY", "NKG7"       # NK cells
)

# Filter markers that exist in the dataset
markers_present <- markers_to_plot[markers_to_plot %in% rownames(pbmc)]

if(length(markers_present) > 0) {
  pdf("known_markers_feature_plot.pdf", width = 15, height = 12)
  print(FeaturePlot(pbmc, features = markers_present, ncol = 3))
  dev.off()
  
  pdf("known_markers_violin_plot.pdf", width = 15, height = 10)
  print(VlnPlot(pbmc, features = markers_present, ncol = 3, pt.size = 0))
  dev.off()
}

# ============================================================================
# STEP 13: Save Final Seurat Object
# ============================================================================
cat("\nSaving final Seurat object...\n")

saveRDS(pbmc, file = "pbmc_final.rds")

cat("\n✓ Analysis complete! Seurat object saved as 'pbmc_final.rds'\n")

# ============================================================================
# STEP 14: Generate Summary Report
# ============================================================================
cat("\n" , rep("=", 70), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Total cells analyzed: %d\n", ncol(pbmc)))
cat(sprintf("Total genes: %d\n", nrow(pbmc)))
cat(sprintf("Number of clusters: %d\n", length(unique(Idents(pbmc)))))
cat(sprintf("PCs used for clustering: %d\n", n_pcs))
cat(sprintf("Clustering resolution: 0.5\n"))
cat("\nCluster sizes:\n")
print(table(Idents(pbmc)))
cat("\nOutput files generated:\n")
cat("  - pbmc_final.rds (Seurat object)\n")
cat("  - qc_plots_before_filtering.pdf\n")
cat("  - qc_scatter_plots.pdf\n")
cat("  - variable_features.pdf\n")
cat("  - pca_plots.pdf\n")
cat("  - elbow_plot.pdf\n")
cat("  - umap_clusters.pdf\n")
cat("  - umap_qc_metrics.pdf\n")
cat("  - cluster_markers.csv\n")
cat("  - top5_markers_per_cluster.csv\n")
cat("  - top_markers_heatmap.pdf\n")
cat("  - known_markers_feature_plot.pdf\n")
cat("  - known_markers_violin_plot.pdf\n")
cat(rep("=", 70), "\n", sep = "")

# ============================================================================
# To reload the object later, use:
# pbmc <- readRDS("pbmc_final.rds")
# ============================================================================
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
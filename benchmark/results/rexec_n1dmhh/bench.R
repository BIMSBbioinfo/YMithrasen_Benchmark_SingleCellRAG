setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Complete Seurat scRNA-seq Analysis Pipeline for PBMC Data
# Production-ready script with comprehensive QC, analysis, and output generation

```r
#!/usr/bin/env Rscript

# =============================================================================
# PBMC scRNA-seq Analysis Pipeline using Seurat
# =============================================================================
# Description: Complete workflow from raw counts to cell type annotation
# Author: BioAssistant (ROLV)
# Date: 2024
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

# Set up parallel processing
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input/Output paths
INPUT_DIR <- "data/raw_counts"
OUTPUT_DIR <- "results/seurat_analysis"
QC_DIR <- file.path(OUTPUT_DIR, "qc_plots")
ANALYSIS_DIR <- file.path(OUTPUT_DIR, "analysis")
MARKERS_DIR <- file.path(OUTPUT_DIR, "markers")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MARKERS_DIR, recursive = TRUE, showWarnings = FALSE)

# QC thresholds
MIN_FEATURES <- 200
MAX_FEATURES <- 2500
MAX_MT_PERCENT <- 5

# Analysis parameters
N_PCS <- 30
RESOLUTION <- 0.5

# =============================================================================
# LOGGING FUNCTION
# =============================================================================

log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message))
}

# =============================================================================
# 1. DATA LOADING
# =============================================================================

log_message("Starting PBMC scRNA-seq analysis pipeline")
log_message("Loading raw count data...")

# Load data (adjust based on your input format)
# Option 1: 10X format
if (dir.exists(file.path(INPUT_DIR, "filtered_feature_bc_matrix"))) {
  pbmc.data <- Read10X(data.dir = file.path(INPUT_DIR, "filtered_feature_bc_matrix"))
  pbmc <- CreateSeuratObject(counts = pbmc.data, 
                             project = "PBMC", 
                             min.cells = 3, 
                             min.features = MIN_FEATURES)
}

# Option 2: Matrix file
# pbmc.data <- read.table(file.path(INPUT_DIR, "counts_matrix.txt"), 
#                         header = TRUE, row.names = 1)
# pbmc <- CreateSeuratObject(counts = pbmc.data, 
#                           project = "PBMC", 
#                           min.cells = 3, 
#                           min.features = MIN_FEATURES)

log_message(sprintf("Loaded %d cells and %d genes", 
                   ncol(pbmc), nrow(pbmc)))

# =============================================================================
# 2. QUALITY CONTROL
# =============================================================================

log_message("Performing quality control...")

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate ribosomal percentage
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# QC metrics before filtering
qc_before <- data.frame(
  nCount_RNA = pbmc$nCount_RNA,
  nFeature_RNA = pbmc$nFeature_RNA,
  percent.mt = pbmc$percent.mt,
  percent.ribo = pbmc$percent.ribo
)

# Violin plots - before filtering
p1 <- VlnPlot(pbmc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
              ncol = 4, pt.size = 0.1) +
  plot_annotation(title = "QC Metrics Before Filtering")

ggsave(file.path(QC_DIR, "qc_violin_before.png"), 
       plot = p1, width = 16, height = 5, dpi = 300)

# Feature scatter plots
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = MAX_MT_PERCENT, linetype = "dashed", color = "red")

p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(MIN_FEATURES, MAX_FEATURES), 
             linetype = "dashed", color = "red")

ggsave(file.path(QC_DIR, "qc_scatter.png"), 
       plot = p2 + p3, width = 12, height = 5, dpi = 300)

# Apply QC filters
log_message(sprintf("Applying QC filters: nFeature %d-%d, MT < %d%%", 
                   MIN_FEATURES, MAX_FEATURES, MAX_MT_PERCENT))

pbmc <- subset(pbmc, 
               subset = nFeature_RNA > MIN_FEATURES & 
                       nFeature_RNA < MAX_FEATURES & 
                       percent.mt < MAX_MT_PERCENT)

log_message(sprintf("After filtering: %d cells retained", ncol(pbmc)))

# QC metrics after filtering
p4 <- VlnPlot(pbmc, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
              ncol = 4, pt.size = 0.1) +
  plot_annotation(title = "QC Metrics After Filtering")

ggsave(file.path(QC_DIR, "qc_violin_after.png"), 
       plot = p4, width = 16, height = 5, dpi = 300)

# Save QC summary
qc_summary <- data.frame(
  Metric = c("Total cells (before)", "Total cells (after)", 
             "Cells removed", "Percent retained"),
  Value = c(nrow(qc_before), ncol(pbmc), 
            nrow(qc_before) - ncol(pbmc),
            round(ncol(pbmc)/nrow(qc_before)*100, 2))
)

write.csv(qc_summary, file.path(QC_DIR, "qc_summary.csv"), row.names = FALSE)

# =============================================================================
# 3. NORMALIZATION
# =============================================================================

log_message("Normalizing data...")

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# =============================================================================
# 4. FEATURE SELECTION
# =============================================================================

log_message("Identifying highly variable features...")

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Plot variable features
top10 <- head(VariableFeatures(pbmc), 10)

p5 <- VariableFeaturePlot(pbmc)
p6 <- LabelPoints(plot = p5, points = top10, repel = TRUE)

ggsave(file.path(ANALYSIS_DIR, "variable_features.png"), 
       plot = p6, width = 10, height = 6, dpi = 300)

# Save variable features
write.csv(data.frame(Gene = VariableFeatures(pbmc)),
          file.path(ANALYSIS_DIR, "variable_features.csv"),
          row.names = FALSE)

# =============================================================================
# 5. SCALING AND PCA
# =============================================================================

log_message("Scaling data and running PCA...")

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = N_PCS)

# Visualize PCA
p7 <- DimPlot(pbmc, reduction = "pca") + 
  ggtitle("PCA Plot")

p8 <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, ncol = 3)

ggsave(file.path(ANALYSIS_DIR, "pca_plot.png"), 
       plot = p7, width = 8, height = 6, dpi = 300)

ggsave(file.path(ANALYSIS_DIR, "pca_heatmap.png"), 
       plot = p8, width = 12, height = 16, dpi = 300)

# Elbow plot for PC selection
p9 <- ElbowPlot(pbmc, ndims = N_PCS) +
  ggtitle("Elbow Plot - PC Selection")

ggsave(file.path(ANALYSIS_DIR, "elbow_plot.png"), 
       plot = p9, width = 8, height = 6, dpi = 300)

# =============================================================================
# 6. CLUSTERING
# =============================================================================

log_message("Performing clustering...")

# Determine optimal number of PCs (using elbow method)
# You can adjust this based on the elbow plot
n_pcs_use <- 20

pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs_use)
pbmc <- FindClusters(pbmc, resolution = RESOLUTION)

log_message(sprintf("Identified %d clusters", 
                   length(unique(Idents(pbmc)))))

# =============================================================================
# 7. UMAP DIMENSIONALITY REDUCTION
# =============================================================================

log_message("Running UMAP...")

pbmc <- RunUMAP(pbmc, dims = 1:n_pcs_use)

# UMAP visualization
p10 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("UMAP - Cell Clusters")

ggsave(file.path(ANALYSIS_DIR, "umap_clusters.png"), 
       plot = p10, width = 10, height = 8, dpi = 300)

# UMAP colored by QC metrics
p11 <- FeaturePlot(pbmc, 
                   features = c("nFeature_RNA", "nCount_RNA", 
                               "percent.mt", "percent.ribo"),
                   ncol = 2)

ggsave(file.path(ANALYSIS_DIR, "umap_qc_metrics.png"), 
       plot = p11, width = 12, height = 10, dpi = 300)

# =============================================================================
# 8. FIND CLUSTER MARKERS
# =============================================================================

log_message("Finding cluster markers...")

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save all markers
write.csv(pbmc.markers, 
          file.path(MARKERS_DIR, "all_cluster_markers.csv"),
          row.names = FALSE)

write.csv(top_markers, 
          file.path(MARKERS_DIR, "top10_markers_per_cluster.csv"),
          row.names = FALSE)

# Heatmap of top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p12 <- DoHeatmap(pbmc, features = top5$gene) + 
  NoLegend() +
  ggtitle("Top 5 Markers per Cluster")

ggsave(file.path(MARKERS_DIR, "marker_heatmap.png"), 
       plot = p12, width = 14, height = 12, dpi = 300)

# =============================================================================
# 9. CELL TYPE ANNOTATION
# =============================================================================

log_message("Annotating cell types...")

# Known PBMC markers
pbmc_markers <- c(
  "CD3D", "CD3E",      # T cells
  "CD4",               # CD4+ T cells
  "CD8A", "CD8B",      # CD8+ T cells
  "IL7R", "CCR7",      # Naive T cells
  "GZMK",              # Memory T cells
  "GNLY", "NKG7",      # NK cells
  "MS4A1", "CD79A",    # B cells
  "CD14", "LYZ",       # Monocytes
  "FCGR3A", "MS4A7",   # FCGR3A+ Monocytes
  "FCER1A", "CST3",    # Dendritic cells
  "PPBP"               # Platelets
)

# Filter markers present in dataset
pbmc_markers_present <- pbmc_markers[pbmc_markers %in% rownames(pbmc)]

# Feature plots for key markers
p13 <- FeaturePlot(pbmc, 
                   features = pbmc_markers_present[1:min(9, length(pbmc_markers_present))],
                   ncol = 3)

ggsave(file.path(ANALYSIS_DIR, "marker_expression_umap.png"), 
       plot = p13, width = 15, height = 12, dpi = 300)

# Dot plot of markers
p14 <- DotPlot(pbmc, features = pbmc_markers_present) + 
  RotatedAxis() +
  ggtitle("Marker Expression by Cluster")

ggsave(file.path(ANALYSIS_DIR, "marker_dotplot.png"), 
       plot = p14, width = 12, height = 6, dpi = 300)

# Violin plots of key markers
p15 <- VlnPlot(pbmc, 
               features = pbmc_markers_present[1:min(6, length(pbmc_markers_present))],
               ncol = 3, pt.size = 0)

ggsave(file.path(ANALYSIS_DIR, "marker_violin.png"), 
       plot = p15, width = 15, height = 8, dpi = 300)

# Manual annotation (adjust based on your marker expression)
new.cluster.ids <- c(
  "0" = "CD4 T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "Dendritic cells",
  "7" = "Platelets"
)

# Apply annotations (only for clusters that exist)
existing_clusters <- levels(Idents(pbmc))
annotations <- new.cluster.ids[existing_clusters]
annotations[is.na(annotations)] <- paste0("Cluster ", existing_clusters[is.na(annotations)])

names(annotations) <- existing_clusters
pbmc <- RenameIdents(pbmc, annotations)

# UMAP with cell type labels
p16 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("UMAP - Annotated Cell Types")

ggsave(file.path(ANALYSIS_DIR, "umap_annotated.png"), 
       plot = p16, width = 12, height = 8, dpi = 300)

# =============================================================================
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
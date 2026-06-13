setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Complete Seurat scRNA-seq Analysis Pipeline for PBMC Data
# Production-ready script with comprehensive QC, analysis, and visualization

```r
#!/usr/bin/env Rscript

# =============================================================================
# PBMC scRNA-seq Analysis Pipeline using Seurat
# Production-ready workflow with quality control and visualization
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set parameters
config <- list(
  # Input/Output
  input_dir = "data/raw_counts",  # Directory containing count matrix
  output_dir = "results/seurat_analysis",
  project_name = "PBMC_scRNA",
  
  # QC thresholds
  min_features = 200,      # Minimum genes per cell
  max_features = 2500,     # Maximum genes per cell
  max_mt_percent = 5,      # Maximum mitochondrial percentage
  
  # Analysis parameters
  n_variable_features = 2000,
  n_pcs = 30,
  resolution = 0.5,
  n_neighbors = 30,
  
  # Visualization
  fig_width = 10,
  fig_height = 8,
  dpi = 300
)

# Create output directories
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(config$output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(config$output_dir, "objects"), showWarnings = FALSE)

# Set up logging
log_file <- file.path(config$output_dir, "analysis_log.txt")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message <- paste0("[", timestamp, "] ", msg)
  cat(message, "\n")
  cat(message, "\n", file = log_file, append = TRUE)
}

log_message("Starting Seurat analysis pipeline")

# =============================================================================
# 1. DATA LOADING
# =============================================================================

log_message("Loading count matrix...")

# Load data (adjust based on your input format)
# Option 1: 10X format
if (dir.exists(file.path(config$input_dir, "filtered_feature_bc_matrix"))) {
  counts <- Read10X(data.dir = file.path(config$input_dir, "filtered_feature_bc_matrix"))
} else if (file.exists(file.path(config$input_dir, "matrix.mtx"))) {
  # Option 2: Matrix Market format
  counts <- readMM(file.path(config$input_dir, "matrix.mtx"))
  genes <- read.table(file.path(config$input_dir, "genes.tsv"), 
                      header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(file.path(config$input_dir, "barcodes.tsv"), 
                         header = FALSE, stringsAsFactors = FALSE)
  rownames(counts) <- genes$V2
  colnames(counts) <- barcodes$V1
} else {
  # Option 3: CSV/TSV format
  counts <- read.csv(file.path(config$input_dir, "counts.csv"), row.names = 1)
}

log_message(paste("Loaded", nrow(counts), "genes and", ncol(counts), "cells"))

# =============================================================================
# 2. CREATE SEURAT OBJECT
# =============================================================================

log_message("Creating Seurat object...")

pbmc <- CreateSeuratObject(
  counts = counts,
  project = config$project_name,
  min.cells = 3,
  min.features = config$min_features
)

log_message(paste("Seurat object created with", ncol(pbmc), "cells"))

# =============================================================================
# 3. QUALITY CONTROL
# =============================================================================

log_message("Performing quality control...")

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate ribosomal percentage
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# QC metrics before filtering
qc_before <- data.frame(
  n_cells = ncol(pbmc),
  median_genes = median(pbmc$nFeature_RNA),
  median_UMIs = median(pbmc$nCount_RNA),
  median_mt = median(pbmc$percent.mt)
)

# Visualize QC metrics
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)

ggsave(
  filename = file.path(config$output_dir, "figures", "01_QC_violin_before.png"),
  plot = p1,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# Feature scatter plots
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

p_scatter <- p2 + p3
ggsave(
  filename = file.path(config$output_dir, "figures", "02_QC_scatter.png"),
  plot = p_scatter,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# Filter cells
log_message("Filtering cells based on QC metrics...")
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > config$min_features & 
           nFeature_RNA < config$max_features & 
           percent.mt < config$max_mt_percent
)

log_message(paste("Retained", ncol(pbmc), "cells after filtering"))

# QC metrics after filtering
qc_after <- data.frame(
  n_cells = ncol(pbmc),
  median_genes = median(pbmc$nFeature_RNA),
  median_UMIs = median(pbmc$nCount_RNA),
  median_mt = median(pbmc$percent.mt)
)

# Save QC summary
qc_summary <- rbind(
  cbind(stage = "before_filtering", qc_before),
  cbind(stage = "after_filtering", qc_after)
)
write.csv(qc_summary, 
          file.path(config$output_dir, "tables", "QC_summary.csv"),
          row.names = FALSE)

# =============================================================================
# 4. NORMALIZATION
# =============================================================================

log_message("Normalizing data...")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# =============================================================================
# 5. FEATURE SELECTION
# =============================================================================

log_message("Identifying highly variable features...")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = config$n_variable_features
)

# Plot variable features
top10 <- head(VariableFeatures(pbmc), 10)
p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE)

ggsave(
  filename = file.path(config$output_dir, "figures", "03_variable_features.png"),
  plot = p5,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# Save variable features
write.csv(
  data.frame(gene = VariableFeatures(pbmc)),
  file.path(config$output_dir, "tables", "variable_features.csv"),
  row.names = FALSE
)

# =============================================================================
# 6. SCALING AND PCA
# =============================================================================

log_message("Scaling data and running PCA...")

# Scale data
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)

# Run PCA
pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = config$n_pcs
)

# Visualize PCA
p6 <- DimPlot(pbmc, reduction = "pca")
ggsave(
  filename = file.path(config$output_dir, "figures", "04_PCA_plot.png"),
  plot = p6,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# PCA loadings
p7 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave(
  filename = file.path(config$output_dir, "figures", "05_PCA_loadings.png"),
  plot = p7,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# Elbow plot
p8 <- ElbowPlot(pbmc, ndims = config$n_pcs)
ggsave(
  filename = file.path(config$output_dir, "figures", "06_elbow_plot.png"),
  plot = p8,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# =============================================================================
# 7. CLUSTERING
# =============================================================================

log_message("Performing clustering...")

# Determine dimensionality (you may want to adjust based on elbow plot)
dims_use <- 1:20

# Find neighbors
pbmc <- FindNeighbors(
  pbmc,
  dims = dims_use,
  k.param = config$n_neighbors
)

# Find clusters
pbmc <- FindClusters(
  pbmc,
  resolution = config$resolution,
  algorithm = 1
)

log_message(paste("Identified", length(unique(Idents(pbmc))), "clusters"))

# =============================================================================
# 8. UMAP VISUALIZATION
# =============================================================================

log_message("Running UMAP...")

pbmc <- RunUMAP(
  pbmc,
  dims = dims_use,
  n.neighbors = config$n_neighbors
)

# UMAP by cluster
p9 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP - Clusters")
ggsave(
  filename = file.path(config$output_dir, "figures", "07_UMAP_clusters.png"),
  plot = p9,
  width = config$fig_width,
  height = config$fig_height,
  dpi = config$dpi
)

# UMAP by QC metrics
p10 <- FeaturePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave(
  filename = file.path(config$output_dir, "figures", "08_UMAP_QC_metrics.png"),
  plot = p10,
  width = config$fig_width * 1.5,
  height = config$fig_height,
  dpi = config$dpi
)

# =============================================================================
# 9. FIND CLUSTER MARKERS
# =============================================================================

log_message("Finding cluster markers...")

# Find markers for all clusters
pbmc_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Save all markers
write.csv(
  pbmc_markers,
  file.path(config$output_dir, "tables", "all_cluster_markers.csv"),
  row.names = FALSE
)

# Get top markers per cluster
top_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(
  top_markers,
  file.path(config$output_dir, "tables", "top10_markers_per_cluster.csv"),
  row.names = FALSE
)

# Heatmap of top markers
top5 <- pbmc_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p11 <- DoHeatmap(pbmc, features = top5$gene) + NoLegend()
ggsave(
  filename = file.path(config$output_dir, "figures", "09_marker_heatmap.png"),
  plot = p11,
  width = config$fig_width * 1.2,
  height = config$fig_height * 1.5,
  dpi = config$dpi
)

# =============================================================================
# 10. CELL TYPE ANNOTATION (PBMC-specific)
# =============================================================================

log_message("Annotating cell types...")

# Known PBMC markers
pbmc_markers_known <- c(
  "CD3D", "CD3E",      # T cells
  "CD4",               # CD4+ T cells
  "CD8A", "CD8B",      # CD8+ T cells
  "MS4A1", "CD79A",    # B cells
  "CD14", "LYZ",       # Monocytes
  "FCGR3A", "MS4A7",   # FCGR3A+ Monocytes
  "GNLY", "NKG7",      # NK cells
  "FCER1A", "CST3"     # Dendritic cells
)

# Feature plots for known markers
p12 <- FeaturePlot(
  pbmc,
  features = pbmc_markers_known[pbmc_markers_known %in% rownames(pbmc)],
  ncol = 4
)
ggsave(
  filename = file.path(config$output_dir, "figures", "10_known_markers.png"),
  plot = p12,
  width = config$fig_width * 2,
  height = config$fig_height * 2,
  dpi = config$dpi
)

# Dot plot of markers
p13 <- DotPlot(
  pbmc,
  features = pbmc_markers_known[pbmc_markers_known %in% rownames(pbmc)]
) + RotatedAxis()
ggsave(
  filename = file.path(config$output_dir, "figures", "11_marker_dotplot.png"),
  plot = p13,
  width = config$fig_width * 1.2,
  height = config$fig_height,
  dpi = config$dpi
)

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================

log_message("Saving Seurat object and metadata...")

# Save Seurat object
saveRDS(
  pbmc,
  file = file.path(config$output_dir, "objects", "pbmc_seurat_object.rds")
)

# Save metadata
metadata <- pbmc@meta.data
metadata$barcode <- rownames(metadata)
write.csv(
  metadata,
  file.path(config$output_dir, "tables", "cell_metadata.csv"),
  row.names = FALSE
)

# Save cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(
  cluster_assignments,
  file.path(config$output_dir, "tables", "cluster_assignments.csv"),
  row.names = FALSE
)

# Save UMAP coordinates
umap_coords <- as.data.frame(Embeddings(pb
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
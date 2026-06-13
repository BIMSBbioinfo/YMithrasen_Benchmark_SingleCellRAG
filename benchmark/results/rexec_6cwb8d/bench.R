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
# PBMC Single-Cell RNA-seq Analysis Pipeline using Seurat
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

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set up parallel processing
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# Define paths
INPUT_DIR <- "data/raw_counts"
OUTPUT_DIR <- "results"
QC_DIR <- file.path(OUTPUT_DIR, "qc")
ANALYSIS_DIR <- file.path(OUTPUT_DIR, "analysis")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ANALYSIS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# QC thresholds
MIN_FEATURES <- 200
MAX_FEATURES <- 2500
MAX_MT_PERCENT <- 5

# Analysis parameters
N_PCS <- 30
RESOLUTION <- 0.5

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

save_plot <- function(plot, filename, width = 10, height = 8) {
  ggsave(
    filename = file.path(FIGURES_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", Sys.time(), message))
}

# =============================================================================
# 1. DATA LOADING
# =============================================================================

log_message("Loading raw count data...")

# Option A: Load from 10X output directory
# pbmc.data <- Read10X(data.dir = INPUT_DIR)

# Option B: Load from matrix file (adjust as needed)
# pbmc.data <- read.table(file.path(INPUT_DIR, "counts.txt"), 
#                         header = TRUE, row.names = 1)

# For demonstration, using example data
# Replace with your actual data loading
pbmc.data <- Read10X(data.dir = INPUT_DIR)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC_scRNA",
  min.cells = 3,
  min.features = MIN_FEATURES
)

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

# Visualize QC metrics before filtering
qc_before <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "QC Metrics Before Filtering")

save_plot(qc_before, "01_qc_before_filtering.png", width = 15, height = 6)

# Feature scatter plots
scatter1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_scatter <- scatter1 + scatter2

save_plot(qc_scatter, "02_qc_scatter.png", width = 12, height = 6)

# Save QC statistics
qc_stats_before <- data.frame(
  Metric = c("Total_Cells", "Mean_Features", "Median_Features", 
             "Mean_Counts", "Median_MT_Percent"),
  Value = c(
    ncol(pbmc),
    mean(pbmc$nFeature_RNA),
    median(pbmc$nFeature_RNA),
    mean(pbmc$nCount_RNA),
    median(pbmc$percent.mt)
  )
)

write.csv(qc_stats_before, 
          file.path(QC_DIR, "qc_stats_before_filtering.csv"),
          row.names = FALSE)

# Filter cells
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > MIN_FEATURES & 
           nFeature_RNA < MAX_FEATURES & 
           percent.mt < MAX_MT_PERCENT
)

log_message(sprintf("After filtering: %d cells retained", ncol(pbmc)))

# QC after filtering
qc_after <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
) + plot_annotation(title = "QC Metrics After Filtering")

save_plot(qc_after, "03_qc_after_filtering.png", width = 15, height = 6)

# =============================================================================
# 3. NORMALIZATION
# =============================================================================

log_message("Normalizing data...")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# =============================================================================
# 4. FEATURE SELECTION
# =============================================================================

log_message("Identifying highly variable features...")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
log_message(paste("Top 10 variable genes:", paste(top10, collapse = ", ")))

# Plot variable features
var_plot <- VariableFeaturePlot(pbmc)
var_plot_labeled <- LabelPoints(
  plot = var_plot,
  points = top10,
  repel = TRUE
)

save_plot(var_plot_labeled, "04_variable_features.png", width = 10, height = 8)

# =============================================================================
# 5. SCALING AND PCA
# =============================================================================

log_message("Scaling data and running PCA...")

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = N_PCS
)

# Visualize PCA results
pca_viz <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
save_plot(pca_viz, "05_pca_loadings.png", width = 10, height = 8)

pca_plot <- DimPlot(pbmc, reduction = "pca")
save_plot(pca_plot, "06_pca_plot.png", width = 8, height = 8)

# Elbow plot for PC selection
elbow <- ElbowPlot(pbmc, ndims = N_PCS)
save_plot(elbow, "07_elbow_plot.png", width = 8, height = 6)

# Heatmap of top PCs
pc_heatmap <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
save_plot(pc_heatmap, "08_pc_heatmap.png", width = 12, height = 15)

# =============================================================================
# 6. CLUSTERING
# =============================================================================

log_message("Performing clustering...")

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:20)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = RESOLUTION)

log_message(sprintf("Identified %d clusters", 
                   length(unique(Idents(pbmc)))))

# =============================================================================
# 7. DIMENSIONALITY REDUCTION (UMAP/tSNE)
# =============================================================================

log_message("Running UMAP and tSNE...")

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Run tSNE
pbmc <- RunTSNE(pbmc, dims = 1:20)

# Visualize clusters
umap_clusters <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP - Clusters")
save_plot(umap_clusters, "09_umap_clusters.png", width = 10, height = 8)

tsne_clusters <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) +
  ggtitle("tSNE - Clusters")
save_plot(tsne_clusters, "10_tsne_clusters.png", width = 10, height = 8)

# =============================================================================
# 8. FIND CLUSTER MARKERS
# =============================================================================

log_message("Finding cluster markers...")

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save all markers
write.csv(
  pbmc.markers,
  file.path(ANALYSIS_DIR, "all_cluster_markers.csv"),
  row.names = FALSE
)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(
  top_markers,
  file.path(ANALYSIS_DIR, "top10_markers_per_cluster.csv"),
  row.names = FALSE
)

# Visualize top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

marker_heatmap <- DoHeatmap(pbmc, features = top5$gene) +
  NoLegend()
save_plot(marker_heatmap, "11_marker_heatmap.png", width = 12, height = 15)

# =============================================================================
# 9. CELL TYPE ANNOTATION
# =============================================================================

log_message("Annotating cell types...")

# Define canonical PBMC markers
canonical_markers <- c(
  "CD3D", "CD3E",      # T cells
  "CD4",               # CD4+ T cells
  "CD8A", "CD8B",      # CD8+ T cells
  "IL7R", "CCR7",      # Naive T cells
  "CD14", "LYZ",       # Monocytes
  "MS4A1", "CD79A",    # B cells
  "GNLY", "NKG7",      # NK cells
  "FCGR3A", "MS4A7",   # FCGR3A+ Monocytes
  "FCER1A", "CST3",    # Dendritic cells
  "PPBP"               # Platelets
)

# Feature plots for canonical markers
for (marker in canonical_markers) {
  if (marker %in% rownames(pbmc)) {
    p <- FeaturePlot(pbmc, features = marker, reduction = "umap")
    save_plot(p, paste0("12_feature_", marker, ".png"), width = 8, height = 7)
  }
}

# Violin plots for key markers
marker_violin <- VlnPlot(
  pbmc,
  features = canonical_markers[canonical_markers %in% rownames(pbmc)],
  stack = TRUE,
  flip = TRUE
) + NoLegend()

save_plot(marker_violin, "13_marker_violin.png", width = 10, height = 12)

# Dot plot for markers
marker_dotplot <- DotPlot(
  pbmc,
  features = canonical_markers[canonical_markers %in% rownames(pbmc)]
) + RotatedAxis()

save_plot(marker_dotplot, "14_marker_dotplot.png", width = 12, height = 6)

# Manual annotation based on markers (adjust based on your data)
new.cluster.ids <- c(
  "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
  "FCGR3A+ Mono", "NK", "DC", "Platelet"
)

# Only rename if number of clusters matches
if (length(unique(Idents(pbmc))) == length(new.cluster.ids)) {
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  pbmc$cell_type <- Idents(pbmc)
  
  # Plot annotated UMAP
  umap_annotated <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
    ggtitle("UMAP - Annotated Cell Types")
  save_plot(umap_annotated, "15_umap_annotated.png", width = 12, height = 8)
}

# =============================================================================
# 10. DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

log_message("Performing differential expression analysis...")

# Example: Compare two cell types (adjust as needed)
if ("cell_type" %in% colnames(pbmc@meta.data)) {
  cell_types <- unique(pbmc$cell_type)
  if (length(cell_types) >= 2) {
    de_results <- FindMarkers(
      pbmc,
      ident.1 = cell_types[1],
      ident.2 = cell_types[2],
      min.pct = 0.25
    )
    
    write.csv(
      de_results,
      file.path(ANALYSIS_DIR, 
                paste0("de_", cell_types[1], "_vs_", cell_types[2], ".csv"))
    )
  }
}

# =============================================================================
# 11. SAVE OUTPUTS
# =============================================================================

log_message("Saving final outputs...")

# Save Seurat object
saveRDS(pbmc, file = file.path(OUTPUT_DIR, "pbmc_seurat_object.rds"))

# Save metadata
write.csv(
  pbmc@meta.data,
  file.path(ANALYSIS_DIR, "cell_metadata.csv")
)

# Save cluster composition
cluster_composition <- table(Idents(pbmc))
write.csv(
  as.data.frame(cluster_composition),
  file.path(ANALYSIS_DIR, "cluster_composition.csv")
)

# Export normalized counts
normalized_counts <- GetAssayData(pbmc, slot = "data")
write.csv(
  as.matrix(normalized_counts),
  file.path(ANALYSIS_DIR, "normalized_counts.csv.gz")
)

# Create summary report
summary_stats <- data.frame(
  Metric = c(
    "Total Cells After QC",
    "Total Genes",
    "Number of Clusters",
    "Number of Variable Features",
    "PCs Used"
  ),
  Value = c(
    ncol(pbmc
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
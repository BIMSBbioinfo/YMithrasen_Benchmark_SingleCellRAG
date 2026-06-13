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
  library(Matrix)
  library(future)
})

# Set up parallel processing
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input/Output paths
INPUT_DIR <- "data/raw_counts"  # Directory with matrix.mtx, genes.tsv, barcodes.tsv
OUTPUT_DIR <- "results/seurat_analysis"
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# QC thresholds
MIN_FEATURES <- 200      # Minimum genes per cell
MAX_FEATURES <- 5000     # Maximum genes per cell
MAX_MT_PERCENT <- 15     # Maximum mitochondrial percentage

# Analysis parameters
N_VARIABLE_FEATURES <- 2000
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

# Load 10X data (adjust based on your input format)
# Option 1: 10X format (matrix.mtx, genes.tsv, barcodes.tsv)
if (file.exists(file.path(INPUT_DIR, "matrix.mtx.gz")) || 
    file.exists(file.path(INPUT_DIR, "matrix.mtx"))) {
  counts <- Read10X(data.dir = INPUT_DIR)
} else if (file.exists(file.path(INPUT_DIR, "counts.rds"))) {
  # Option 2: RDS file
  counts <- readRDS(file.path(INPUT_DIR, "counts.rds"))
} else {
  stop("No valid input data found. Please provide matrix.mtx or counts.rds")
}

log_message(sprintf("Loaded data: %d genes x %d cells", 
                   nrow(counts), ncol(counts)))

# =============================================================================
# 2. CREATE SEURAT OBJECT
# =============================================================================

log_message("Creating Seurat object...")

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "PBMC_scRNAseq",
  min.cells = 3,
  min.features = MIN_FEATURES
)

log_message(sprintf("Seurat object created: %d cells retained", ncol(pbmc)))

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
  nCells = ncol(pbmc),
  median_nFeature = median(pbmc$nFeature_RNA),
  median_nCount = median(pbmc$nCount_RNA),
  median_percent_mt = median(pbmc$percent.mt)
)

# Visualize QC metrics
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
ggsave(file.path(FIGURES_DIR, "01_QC_violin_before_filtering.png"), 
       p1, width = 12, height = 4, dpi = 300)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = MAX_MT_PERCENT, linetype = "dashed", color = "red")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(MIN_FEATURES, MAX_FEATURES), 
             linetype = "dashed", color = "red")
p_scatter <- p2 | p3
ggsave(file.path(FIGURES_DIR, "02_QC_scatter_before_filtering.png"), 
       p_scatter, width = 12, height = 5, dpi = 300)

# Apply QC filters
pbmc <- subset(pbmc, subset = 
                 nFeature_RNA > MIN_FEATURES & 
                 nFeature_RNA < MAX_FEATURES & 
                 percent.mt < MAX_MT_PERCENT)

log_message(sprintf("After QC filtering: %d cells retained", ncol(pbmc)))

# QC metrics after filtering
qc_after <- data.frame(
  nCells = ncol(pbmc),
  median_nFeature = median(pbmc$nFeature_RNA),
  median_nCount = median(pbmc$nCount_RNA),
  median_percent_mt = median(pbmc$percent.mt)
)

# Save QC summary
qc_summary <- rbind(
  cbind(Stage = "Before_Filtering", qc_before),
  cbind(Stage = "After_Filtering", qc_after)
)
write.csv(qc_summary, file.path(DATA_DIR, "QC_summary.csv"), row.names = FALSE)

# =============================================================================
# 4. NORMALIZATION
# =============================================================================

log_message("Normalizing data...")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

# =============================================================================
# 5. FEATURE SELECTION
# =============================================================================

log_message("Identifying highly variable features...")

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = N_VARIABLE_FEATURES)

# Plot variable features
top10 <- head(VariableFeatures(pbmc), 10)
p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(file.path(FIGURES_DIR, "03_variable_features.png"), 
       p5, width = 10, height = 6, dpi = 300)

# Save variable features
write.csv(data.frame(Gene = VariableFeatures(pbmc)), 
          file.path(DATA_DIR, "variable_features.csv"), row.names = FALSE)

# =============================================================================
# 6. SCALING AND PCA
# =============================================================================

log_message("Scaling data and running PCA...")

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), 
              npcs = N_PCS, verbose = FALSE)

# Visualize PCA
p6 <- DimPlot(pbmc, reduction = "pca") + ggtitle("PCA")
ggsave(file.path(FIGURES_DIR, "04_PCA_plot.png"), 
       p6, width = 8, height = 6, dpi = 300)

p7 <- ElbowPlot(pbmc, ndims = N_PCS)
ggsave(file.path(FIGURES_DIR, "05_PCA_elbow_plot.png"), 
       p7, width = 8, height = 6, dpi = 300)

# PCA loadings
p8 <- VizDimLoadings(pbmc, dims = 1:4, reduction = "pca", ncol = 2)
ggsave(file.path(FIGURES_DIR, "06_PCA_loadings.png"), 
       p8, width = 12, height = 10, dpi = 300)

# =============================================================================
# 7. CLUSTERING
# =============================================================================

log_message("Performing clustering...")

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = RESOLUTION)

log_message(sprintf("Identified %d clusters", 
                   length(unique(Idents(pbmc)))))

# =============================================================================
# 8. UMAP DIMENSIONALITY REDUCTION
# =============================================================================

log_message("Running UMAP...")

pbmc <- RunUMAP(pbmc, dims = 1:20)

# UMAP visualization
p9 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) + 
  ggtitle("UMAP - Cell Clusters")
ggsave(file.path(FIGURES_DIR, "07_UMAP_clusters.png"), 
       p9, width = 10, height = 8, dpi = 300)

# UMAP with QC metrics
p10 <- FeaturePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   ncol = 3)
ggsave(file.path(FIGURES_DIR, "08_UMAP_QC_metrics.png"), 
       p10, width = 15, height = 5, dpi = 300)

# =============================================================================
# 9. FIND CLUSTER MARKERS
# =============================================================================

log_message("Finding cluster markers...")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(pbmc.markers, file.path(DATA_DIR, "all_cluster_markers.csv"), 
          row.names = FALSE)
write.csv(top_markers, file.path(DATA_DIR, "top10_cluster_markers.csv"), 
          row.names = FALSE)

# Heatmap of top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p11 <- DoHeatmap(pbmc, features = top5$gene) + 
  theme(axis.text.y = element_text(size = 6))
ggsave(file.path(FIGURES_DIR, "09_marker_heatmap.png"), 
       p11, width = 12, height = 14, dpi = 300)

# =============================================================================
# 10. CELL TYPE ANNOTATION
# =============================================================================

log_message("Annotating cell types...")

# PBMC canonical markers
canonical_markers <- list(
  "CD14+ Monocytes" = c("CD14", "LYZ"),
  "FCGR3A+ Monocytes" = c("FCGR3A", "MS4A7"),
  "CD4+ T cells" = c("CD4", "IL7R"),
  "CD8+ T cells" = c("CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7"),
  "B cells" = c("MS4A1", "CD79A"),
  "Dendritic cells" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP")
)

# Plot canonical markers
for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  markers <- markers[markers %in% rownames(pbmc)]
  
  if (length(markers) > 0) {
    p <- FeaturePlot(pbmc, features = markers, ncol = length(markers))
    filename <- gsub("[+ ]", "_", cell_type)
    ggsave(file.path(FIGURES_DIR, sprintf("10_markers_%s.png", filename)), 
           p, width = 4 * length(markers), height = 4, dpi = 300)
  }
}

# Dot plot of canonical markers
all_markers <- unlist(canonical_markers)
all_markers <- unique(all_markers[all_markers %in% rownames(pbmc)])

p12 <- DotPlot(pbmc, features = all_markers) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIGURES_DIR, "11_canonical_markers_dotplot.png"), 
       p12, width = 12, height = 6, dpi = 300)

# Manual annotation based on markers (adjust based on your data)
new.cluster.ids <- c(
  "0" = "CD4+ T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8+ T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "Dendritic cells",
  "7" = "Platelets"
)

# Apply annotations (adjust cluster numbers based on your results)
current_clusters <- levels(Idents(pbmc))
annotations <- new.cluster.ids[current_clusters]
annotations[is.na(annotations)] <- paste0("Cluster_", current_clusters[is.na(annotations)])
names(annotations) <- current_clusters

pbmc <- RenameIdents(pbmc, annotations)
pbmc$cell_type <- Idents(pbmc)

# Plot annotated UMAP
p13 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) + 
  ggtitle("UMAP - Annotated Cell Types")
ggsave(file.path(FIGURES_DIR, "12_UMAP_annotated.png"), 
       p13, width = 12, height = 8, dpi = 300)

# =============================================================================
# 11. DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# =============================================================================

log_message("Performing differential expression analysis...")

# Find markers for each cell type
celltype.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25)

write.csv(celltype.markers, 
          file.path(DATA_DIR, "celltype_markers.csv"), 
          row.names = FALSE)

# =============================================================================
# 12. SAVE RESULTS
# =============================================================================

log_message("Saving final results...")

# Save Seurat object
saveRDS(pbmc, file.path(DATA_DIR, "pbmc_seurat_final.rds"))

# Save metadata
metadata <- pbmc@meta.data
metadata$barcode <- rownames(metadata)
write.csv(metadata, file.path(DATA_DIR, "cell_metadata.csv"), row.names
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
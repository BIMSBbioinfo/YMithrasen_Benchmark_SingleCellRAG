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
  message(sprintf("Saved: %s", filename))
}

log_step <- function(message) {
  cat(sprintf("\n[%s] %s\n", Sys.time(), message))
}

# =============================================================================
# 1. DATA LOADING
# =============================================================================

log_step("Loading raw count data...")

# Load 10X data (adjust based on your input format)
# Option 1: 10X format (matrix.mtx, genes.tsv, barcodes.tsv)
if (file.exists(file.path(INPUT_DIR, "matrix.mtx.gz")) || 
    file.exists(file.path(INPUT_DIR, "matrix.mtx"))) {
  counts <- Read10X(data.dir = INPUT_DIR)
} else if (file.exists(file.path(INPUT_DIR, "counts.rds"))) {
  # Option 2: R object
  counts <- readRDS(file.path(INPUT_DIR, "counts.rds"))
} else {
  stop("No valid input data found. Please check INPUT_DIR.")
}

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = counts,
  project = "PBMC_scRNA",
  min.cells = 3,
  min.features = MIN_FEATURES
)

log_step(sprintf("Loaded %d cells and %d genes", ncol(pbmc), nrow(pbmc)))

# =============================================================================
# 2. QUALITY CONTROL
# =============================================================================

log_step("Performing quality control...")

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate ribosomal percentage
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# QC metrics before filtering
qc_before <- data.frame(
  nCount_RNA = pbmc$nCount_RNA,
  nFeature_RNA = pbmc$nFeature_RNA,
  percent.mt = pbmc$percent.mt,
  status = "Before"
)

# Visualize QC metrics
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.1)
save_plot(p1, "01_qc_violin_before_filtering.png", width = 15, height = 6)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = MAX_MT_PERCENT, linetype = "dashed", color = "red")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = c(MIN_FEATURES, MAX_FEATURES), 
             linetype = "dashed", color = "red")
save_plot(p2 + p3, "02_qc_scatter_before_filtering.png", width = 14, height = 6)

# Filter cells
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > MIN_FEATURES & 
           nFeature_RNA < MAX_FEATURES & 
           percent.mt < MAX_MT_PERCENT
)

log_step(sprintf("After filtering: %d cells retained", ncol(pbmc)))

# QC metrics after filtering
qc_after <- data.frame(
  nCount_RNA = pbmc$nCount_RNA,
  nFeature_RNA = pbmc$nFeature_RNA,
  percent.mt = pbmc$percent.mt,
  status = "After"
)

# Save QC summary
qc_summary <- rbind(
  data.frame(
    Metric = c("Total Cells", "Median Genes", "Median UMIs", "Median MT%"),
    Before = c(
      nrow(qc_before),
      median(qc_before$nFeature_RNA),
      median(qc_before$nCount_RNA),
      round(median(qc_before$percent.mt), 2)
    ),
    After = c(
      nrow(qc_after),
      median(qc_after$nFeature_RNA),
      median(qc_after$nCount_RNA),
      round(median(qc_after$percent.mt), 2)
    )
  )
)
write.csv(qc_summary, file.path(DATA_DIR, "qc_summary.csv"), row.names = FALSE)

# =============================================================================
# 3. NORMALIZATION
# =============================================================================

log_step("Normalizing data...")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# =============================================================================
# 4. FEATURE SELECTION
# =============================================================================

log_step("Identifying highly variable features...")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = N_VARIABLE_FEATURES
)

# Plot variable features
top10 <- head(VariableFeatures(pbmc), 10)
p4 <- VariableFeaturePlot(pbmc)
p5 <- LabelPoints(plot = p4, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
save_plot(p5, "03_variable_features.png", width = 12, height = 8)

# =============================================================================
# 5. SCALING AND PCA
# =============================================================================

log_step("Scaling data and running PCA...")

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = N_PCS,
  verbose = FALSE
)

# Visualize PCA
p6 <- DimPlot(pbmc, reduction = "pca") + ggtitle("PCA")
save_plot(p6, "04_pca_plot.png")

p7 <- ElbowPlot(pbmc, ndims = N_PCS)
save_plot(p7, "05_elbow_plot.png")

# PCA heatmaps
pdf(file.path(FIGURES_DIR, "06_pca_heatmaps.pdf"), width = 12, height = 20)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, ncol = 3)
dev.off()

# =============================================================================
# 6. CLUSTERING
# =============================================================================

log_step("Performing clustering...")

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:20)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = RESOLUTION)

log_step(sprintf("Identified %d clusters", length(unique(pbmc$seurat_clusters))))

# =============================================================================
# 7. DIMENSIONALITY REDUCTION (UMAP & tSNE)
# =============================================================================

log_step("Running UMAP and tSNE...")

pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)

# Visualize clusters
p8 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP - Clusters")
p9 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
  ggtitle("tSNE - Clusters")
save_plot(p8 + p9, "07_umap_tsne_clusters.png", width = 16, height = 7)

# =============================================================================
# 8. FIND CLUSTER MARKERS
# =============================================================================

log_step("Finding cluster markers...")

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Save all markers
write.csv(
  pbmc.markers,
  file.path(DATA_DIR, "all_cluster_markers.csv"),
  row.names = FALSE
)

# Top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(
  top_markers,
  file.path(DATA_DIR, "top10_markers_per_cluster.csv"),
  row.names = FALSE
)

# Heatmap of top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

p10 <- DoHeatmap(pbmc, features = top5$gene) + 
  theme(axis.text.y = element_text(size = 6))
save_plot(p10, "08_top_markers_heatmap.png", width = 14, height = 12)

# =============================================================================
# 9. CELL TYPE ANNOTATION
# =============================================================================

log_step("Annotating cell types...")

# PBMC canonical markers
canonical_markers <- list(
  "CD14+ Monocytes" = c("CD14", "LYZ"),
  "FCGR3A+ Monocytes" = c("FCGR3A", "MS4A7"),
  "CD4+ T cells" = c("CD4", "IL7R", "CCR7"),
  "CD8+ T cells" = c("CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7"),
  "B cells" = c("MS4A1", "CD79A"),
  "DC" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP")
)

# Visualize canonical markers
for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  markers <- markers[markers %in% rownames(pbmc)]
  
  if (length(markers) > 0) {
    p <- FeaturePlot(pbmc, features = markers, ncol = length(markers))
    filename <- sprintf("09_markers_%s.png", gsub("[+ ]", "_", cell_type))
    save_plot(p, filename, width = 4 * length(markers), height = 4)
  }
}

# Dot plot of canonical markers
all_canonical <- unlist(canonical_markers)
all_canonical <- all_canonical[all_canonical %in% rownames(pbmc)]

p11 <- DotPlot(pbmc, features = unique(all_canonical)) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p11, "10_canonical_markers_dotplot.png", width = 14, height = 8)

# Manual annotation based on markers (adjust based on your data)
new.cluster.ids <- c(
  "0" = "CD4+ T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8+ T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "DC",
  "7" = "Platelets"
)

# Only rename clusters that exist
existing_clusters <- levels(pbmc$seurat_clusters)
new.cluster.ids <- new.cluster.ids[names(new.cluster.ids) %in% existing_clusters]

pbmc$cell_type <- plyr::mapvalues(
  x = pbmc$seurat_clusters,
  from = names(new.cluster.ids),
  to = new.cluster.ids
)

# Visualize annotated cell types
p12 <- DimPlot(pbmc, reduction = "umap", group.by = "cell_type", 
               label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP - Cell Types")
save_plot(p12, "11_umap_cell_types.png", width = 12, height = 8)

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

log_step("Saving final results...")

# Save Seurat object
saveRDS(pbmc, file.path(DATA_DIR, "pbmc_seurat_final.rds"))

# Save metadata
write.csv(
  pbmc@meta.data,
  file.path(DATA_DIR, "cell_metadata.csv"),
  row.names = TRUE
)

# Save normalized counts
write.csv(
  as.matrix(GetAssayData(pbmc, slot = "data")),
  file.path(DATA_DIR, "normalized_counts.csv.gz")
)

# Save UMAP coordinates
umap_coords <- as.data.frame(Embeddings(pbmc, reduction = "umap"))
umap_coords$cell_type <- pbmc$cell_type
umap_coords$cluster <- pbmc$seurat_clusters
write.csv(
  umap_coords,
  file.path(DATA_DIR, "umap_coordinates.csv"),
  row.names = TRUE
)

# Cell type composition
composition <- table(pbmc$cell_type)
composition_df <- data.frame(
  Cell_Type = names(composition),
  Count = as.numeric(composition),
  Percentage = round(as.numeric(composition) / sum(composition) * 100, 2)
)
write.csv(
  composition_df,
  file.path(DATA_DIR, "cell_type_composition.csv"),
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
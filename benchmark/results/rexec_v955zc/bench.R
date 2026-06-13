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
library(ggplot2)
library(patchwork)

# Assuming you have a Seurat object named 'pbmc' with clusters already defined
# If starting from scratch, you'd need preprocessing steps (normalization, scaling, PCA, clustering)

# ============================================
# 1. FIND MARKER GENES FOR EACH CLUSTER
# ============================================

# Find markers for all clusters (one vs all approach)
all_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,           # Only positive markers
  min.pct = 0.25,            # Expressed in at least 25% of cells
  logfc.threshold = 0.25,    # Log fold-change threshold
  test.use = "wilcox"        # Wilcoxon rank sum test
)

# View top markers per cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top10_markers)

# Save marker genes to file
write.csv(all_markers, "pbmc_all_markers.csv", row.names = FALSE)
write.csv(top10_markers, "pbmc_top10_markers.csv", row.names = FALSE)

# ============================================
# 2. VISUALIZE MARKER GENES
# ============================================

# Heatmap of top markers
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5_markers$gene) + 
  NoLegend() +
  theme(axis.text.y = element_text(size = 8))

# Dot plot of canonical markers
canonical_markers <- c(
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

DotPlot(pbmc, features = canonical_markers) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Feature plots for key markers
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7"), 
            ncol = 2)

# Violin plots
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7"), 
        ncol = 2)

# ============================================
# 3. CELL TYPE ANNOTATION
# ============================================

# Based on canonical markers, annotate clusters
# This is an example - adjust based on your actual cluster markers

new_cluster_ids <- c(
  "0" = "CD4+ T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8+ T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "FCGR3A+ Monocytes",
  "6" = "Dendritic cells",
  "7" = "Platelets"
)

# Rename clusters
pbmc$cell_type <- plyr::mapvalues(
  x = pbmc$seurat_clusters,
  from = names(new_cluster_ids),
  to = new_cluster_ids
)

# Set identity to cell types
Idents(pbmc) <- "cell_type"

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# ============================================
# 4. AUTOMATED ANNOTATION (OPTIONAL)
# ============================================

# Using SingleR for automated annotation
library(SingleR)
library(celldex)

# Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Convert Seurat to SingleCellExperiment
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(pbmc)

# Run SingleR
predictions <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# Add predictions to Seurat object
pbmc$singler_labels <- predictions$labels

# Compare manual and automated annotations
table(pbmc$cell_type, pbmc$singler_labels)

# ============================================
# 5. DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# ============================================

# Find markers between specific cell types
cd4_vs_cd8 <- FindMarkers(
  pbmc,
  ident.1 = "CD4+ T cells",
  ident.2 = "CD8+ T cells",
  min.pct = 0.25
)

head(cd4_vs_cd8, 20)

# ============================================
# 6. QUALITY CONTROL METRICS BY CELL TYPE
# ============================================

# Visualize QC metrics by cell type
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0)

# ============================================
# 7. SAVE ANNOTATED OBJECT
# ============================================

saveRDS(pbmc, "pbmc_annotated.rds")

# Export cell type assignments
cell_annotations <- data.frame(
  cell_id = colnames(pbmc),
  cluster = pbmc$seurat_clusters,
  cell_type = pbmc$cell_type,
  singler_annotation = pbmc$singler_labels
)

write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
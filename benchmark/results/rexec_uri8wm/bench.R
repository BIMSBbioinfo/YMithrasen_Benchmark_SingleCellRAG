setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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

# Assuming you have a Seurat object named 'pbmc' with clusters already defined
# If starting from scratch, you'd need preprocessing steps (normalization, scaling, PCA, clustering)

# ============================================
# 1. Find Marker Genes for Each Cluster
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
write.csv(all_markers, "cluster_markers.csv", row.names = FALSE)

# ============================================
# 2. Visualize Top Markers
# ============================================

# Heatmap of top markers
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5_markers$gene) + 
  theme(axis.text.y = element_text(size = 8))

# Dot plot of top markers
DotPlot(pbmc, features = unique(top5_markers$gene)) + 
  RotatedAxis()

# ============================================
# 3. Cell Type Annotation Based on Known Markers
# ============================================

# Define canonical PBMC markers
canonical_markers <- list(
  "CD4 T cells" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8 T cells" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7", "NCAM1", "KLRD1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "DC" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP", "PF4")
)

# Visualize canonical markers
for (cell_type in names(canonical_markers)) {
  print(FeaturePlot(pbmc, features = canonical_markers[[cell_type]], 
                    ncol = 2) + 
          plot_annotation(title = cell_type))
}

# Dot plot of canonical markers
DotPlot(pbmc, features = canonical_markers) + 
  RotatedAxis() +
  ggtitle("Canonical PBMC Markers")

# ============================================
# 4. Automated Cell Type Annotation
# ============================================

# Manual annotation based on top markers
# Examine top markers for each cluster
for (i in unique(all_markers$cluster)) {
  cat("\n=== Cluster", i, "===\n")
  cluster_markers <- all_markers %>% 
    filter(cluster == i) %>% 
    top_n(n = 10, wt = avg_log2FC)
  print(cluster_markers$gene)
}

# Create annotation mapping (example - adjust based on your data)
cluster_annotations <- c(
  "0" = "CD4 T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "CD16+ Monocytes",
  "6" = "Dendritic cells",
  "7" = "Platelets"
)

# Add cell type annotations to Seurat object
pbmc$cell_type <- cluster_annotations[as.character(Idents(pbmc))]

# Visualize annotated cell types
DimPlot(pbmc, reduction = "umap", group.by = "cell_type", label = TRUE, 
        label.size = 4, repel = TRUE) +
  ggtitle("Annotated PBMC Cell Types")

# ============================================
# 5. Alternative: SingleR for Automated Annotation
# ============================================

library(SingleR)
library(celldex)

# Load reference dataset
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

# Visualize SingleR predictions
DimPlot(pbmc, reduction = "umap", group.by = "singler_labels", 
        label = TRUE, repel = TRUE) +
  ggtitle("SingleR Cell Type Predictions")

# ============================================
# 6. Validation and Quality Control
# ============================================

# Check marker expression in annotated clusters
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7"), 
        group.by = "cell_type", ncol = 2)

# Feature plots for key markers
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7"), 
            ncol = 2)

# Summary statistics
table(pbmc$cell_type)
prop.table(table(pbmc$cell_type))

# ============================================
# 7. Export Results
# ============================================

# Save annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")

# Export cell type assignments
cell_annotations <- data.frame(
  cell_id = colnames(pbmc),
  cluster = Idents(pbmc),
  cell_type = pbmc$cell_type,
  singler_prediction = pbmc$singler_labels
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
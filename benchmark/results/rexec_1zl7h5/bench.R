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
# If starting from raw data, you'd need preprocessing steps first

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
# 2. Visualize Top Marker Genes
# ============================================

# Heatmap of top markers
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5_markers$gene) + 
  theme(axis.text.y = element_text(size = 8))

# Feature plots for specific markers
FeaturePlot(pbmc, features = c("CD3D", "CD8A", "CD4", "MS4A1"))

# Violin plots
VlnPlot(pbmc, features = c("CD3D", "CD8A", "CD4", "MS4A1"))

# Dot plot
DotPlot(pbmc, features = unique(top5_markers$gene)) + 
  RotatedAxis()

# ============================================
# 3. Cell Type Annotation Based on Known Markers
# ============================================

# Define canonical PBMC markers
pbmc_markers <- list(
  "CD4+ T cells" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8+ T cells" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7", "NCAM1", "KLRD1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "DC" = c("FCER1A", "CST3", "CD1C"),
  "Platelets" = c("PPBP", "PF4")
)

# Visualize canonical markers
for (cell_type in names(pbmc_markers)) {
  print(FeaturePlot(pbmc, features = pbmc_markers[[cell_type]], 
                    ncol = 2) + 
          plot_annotation(title = cell_type))
}

# ============================================
# 4. Automated Cell Type Annotation
# ============================================

# Manual annotation based on top markers
# Examine top markers for each cluster
for (i in unique(all_markers$cluster)) {
  cat("\n=== Cluster", i, "===\n")
  cluster_markers <- all_markers %>% 
    filter(cluster == i) %>% 
    arrange(desc(avg_log2FC)) %>%
    head(10)
  print(cluster_markers[, c("gene", "avg_log2FC", "pct.1", "pct.2")])
}

# Create annotation based on marker expression
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
pbmc$cell_type <- new_cluster_ids[as.character(Idents(pbmc))]
Idents(pbmc) <- "cell_type"

# Visualize annotated cell types
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# ============================================
# 5. Alternative: Use SingleR for Automated Annotation
# ============================================

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

# Visualize SingleR predictions
DimPlot(pbmc, reduction = "umap", group.by = "singler_labels", 
        label = TRUE, pt.size = 0.5)

# ============================================
# 6. Quality Control for Annotations
# ============================================

# Check marker expression in annotated clusters
DotPlot(pbmc, features = unlist(pbmc_markers), group.by = "cell_type") + 
  RotatedAxis()

# Calculate percentage of cells per cell type
cell_type_counts <- table(pbmc$cell_type)
cell_type_percentages <- prop.table(cell_type_counts) * 100
print(cell_type_percentages)

# ============================================
# 7. Find Markers Between Specific Cell Types
# ============================================

# Example: CD4+ vs CD8+ T cells
cd4_vs_cd8_markers <- FindMarkers(
  pbmc,
  ident.1 = "CD4+ T cells",
  ident.2 = "CD8+ T cells",
  min.pct = 0.25
)

head(cd4_vs_cd8_markers, 20)

# ============================================
# 8. Export Results
# ============================================

# Save annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")

# Export cell type assignments
cell_annotations <- data.frame(
  cell_id = colnames(pbmc),
  cluster = pbmc$seurat_clusters,
  cell_type = pbmc$cell_type,
  singler_prediction = pbmc$singler_labels
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
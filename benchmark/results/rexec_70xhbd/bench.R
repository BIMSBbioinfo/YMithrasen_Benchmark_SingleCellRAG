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

# Assuming you have a Seurat object named 'pbmc' with clusters already defined
# If starting from scratch, you'd need preprocessing steps (normalization, scaling, PCA, clustering)

# ============================================
# 1. Find Marker Genes for Each Cluster
# ============================================

# Find markers for all clusters (one vs all approach)
pbmc.markers <- FindAllMarkers(
  pbmc, 
  only.pos = TRUE,           # Only positive markers
  min.pct = 0.25,            # Gene expressed in at least 25% of cells
  logfc.threshold = 0.25,    # Minimum log fold change
  test.use = "wilcox"        # Wilcoxon rank sum test
)

# View top markers per cluster
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top10_markers)

# Save marker genes to file
write.csv(pbmc.markers, "pbmc_cluster_markers.csv", row.names = FALSE)

# ============================================
# 2. Visualize Top Marker Genes
# ============================================

# Heatmap of top markers
top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5_markers$gene) + 
  theme(axis.text.y = element_text(size = 8))

# Feature plots for specific markers
FeaturePlot(pbmc, features = c("CD3D", "CD8A", "CD4", "MS4A1"))

# Violin plots
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "FCGR3A"))

# Dot plot
DotPlot(pbmc, features = unique(top5_markers$gene)) + 
  RotatedAxis()

# ============================================
# 3. Cell Type Annotation Based on Markers
# ============================================

# Define canonical PBMC markers
pbmc_markers <- list(
  "CD4+ T cells" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8+ T cells" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7", "NCAM1", "KLRD1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "CD14+ Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "FCGR3A+ Monocytes" = c("FCGR3A", "MS4A7", "CD16"),
  "Dendritic cells" = c("FCER1A", "CST3", "CD1C"),
  "Platelets" = c("PPBP", "PF4")
)

# Automated annotation using marker expression
# Calculate average expression of marker sets per cluster
cluster_avg <- AverageExpression(pbmc, return.seurat = FALSE)$RNA

# Function to score each cluster for each cell type
score_cell_types <- function(cluster_data, marker_list) {
  scores <- sapply(marker_list, function(markers) {
    available_markers <- intersect(markers, rownames(cluster_data))
    if(length(available_markers) > 0) {
      mean(cluster_data[available_markers, ], na.rm = TRUE)
    } else {
      0
    }
  })
  return(scores)
}

# Score each cluster
cluster_scores <- apply(cluster_avg, 2, function(x) {
  score_cell_types(matrix(x, ncol = 1, dimnames = list(names(x), NULL)), pbmc_markers)
})

# Assign cell types based on highest score
predicted_types <- colnames(cluster_scores)[apply(cluster_scores, 2, which.max)]
names(predicted_types) <- colnames(cluster_avg)

print("Predicted cell types:")
print(predicted_types)

# ============================================
# 4. Manual Annotation (Recommended)
# ============================================

# Based on marker gene expression, manually assign cell types
# Example annotation (adjust based on your actual cluster markers):
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

# Rename clusters
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$cell_type <- Idents(pbmc)

# ============================================
# 5. Visualization with Cell Type Labels
# ============================================

# UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Cell Types")

# Save annotated object
saveRDS(pbmc, "pbmc_annotated.rds")

# ============================================
# 6. Validate Annotations
# ============================================

# Check expression of known markers in annotated cell types
markers_to_check <- c("CD3D", "CD4", "CD8A", "MS4A1", "CD14", "FCGR3A", "NKG7")

DotPlot(pbmc, features = markers_to_check) + 
  RotatedAxis() +
  ggtitle("Marker Expression by Cell Type")

# Statistical summary
table(pbmc$cell_type)

# Export cell type assignments
cell_annotations <- data.frame(
  cell_id = colnames(pbmc),
  cluster = pbmc$seurat_clusters,
  cell_type = pbmc$cell_type
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
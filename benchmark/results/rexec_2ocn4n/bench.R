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

# Assuming you have a Seurat object named 'pbmc'
# If starting from count matrix:
# pbmc <- CreateSeuratObject(counts = count_matrix, project = "PBMC")

# Standard preprocessing (if not already done)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:15)

# ============================================
# FIND MARKER GENES FOR EACH CLUSTER
# ============================================

# Method 1: Find markers for all clusters (one vs all)
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

# Method 2: Find markers between specific clusters
# Example: Compare cluster 0 vs cluster 1
cluster0_vs_1 <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1)

# ============================================
# VISUALIZE MARKER GENES
# ============================================

# Heatmap of top markers
top5_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(pbmc, features = top5_markers$gene) + 
  theme(axis.text.y = element_text(size = 8))

# Feature plots for specific markers
FeaturePlot(pbmc, features = c("CD3D", "CD8A", "CD4", "CD14", 
                                "FCGR3A", "MS4A1", "GNLY", "NKG7"))

# Violin plots
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY"))

# Dot plot
DotPlot(pbmc, features = c("CD3D", "CD8A", "CD4", "CD14", 
                           "FCGR3A", "MS4A1", "GNLY", "NKG7")) + 
  RotatedAxis()

# ============================================
# CELL TYPE ANNOTATION
# ============================================

# Define canonical PBMC markers
pbmc_markers <- list(
  "CD4+ T cells" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8+ T cells" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7", "NCAM1", "KLRD1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "CD14+ Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "FCGR3A+ Monocytes" = c("FCGR3A", "MS4A7", "CD16"),
  "Dendritic cells" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP", "PF4")
)

# Visualize canonical markers
for (cell_type in names(pbmc_markers)) {
  print(FeaturePlot(pbmc, features = pbmc_markers[[cell_type]], 
                    ncol = 2) + 
          plot_annotation(title = cell_type))
}

# Manual annotation based on marker expression
new_cluster_ids <- c(
  "0" = "CD14+ Monocytes",
  "1" = "CD4+ T cells",
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

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Annotated PBMC Cell Types")

# ============================================
# AUTOMATED ANNOTATION (OPTIONAL)
# ============================================

# Using SingleR for automated annotation
library(SingleR)
library(celldex)

# Get reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()

# Extract counts and run SingleR
pbmc_counts <- GetAssayData(pbmc, slot = "data")
singler_results <- SingleR(
  test = pbmc_counts,
  ref = ref,
  labels = ref$label.main
)

# Add SingleR annotations to Seurat object
pbmc$singler_labels <- singler_results$labels

# Compare manual vs automated annotation
table(pbmc$cell_type, pbmc$singler_labels)

# ============================================
# STATISTICAL SUMMARY
# ============================================

# Cell counts per cluster
table(Idents(pbmc))

# Percentage of cells per cluster
prop.table(table(Idents(pbmc))) * 100

# Average expression of key markers per cluster
avg_expression <- AverageExpression(pbmc, 
                                   features = c("CD3D", "CD14", "MS4A1", "GNLY"),
                                   return.seurat = FALSE)
print(avg_expression$RNA)

# ============================================
# EXPORT RESULTS
# ============================================

# Save annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")

# Export cell type assignments
cell_annotations <- data.frame(
  cell_id = colnames(pbmc),
  cluster = pbmc$seurat_clusters,
  cell_type = pbmc$cell_type,
  umap1 = Embeddings(pbmc, "umap")[,1],
  umap2 = Embeddings(pbmc, "umap")[,2]
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)

# Generate summary report
cat("=== PBMC scRNA-seq Analysis Summary ===\n")
cat("Total cells:", ncol(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("\nCell type distribution:\n")
print(table(Idents(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
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
# If not, here's a basic preprocessing workflow:

# Load data (example)
# pbmc <- Read10X(data.dir = ".")
# pbmc <- CreateSeuratObject(counts = pbmc, project = "PBMC", min.cells = 3, min.features = 200)

# Standard preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# ============================================
# MARKER GENE IDENTIFICATION
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

# Method 2: Find markers between specific clusters
# cluster1_vs_2 <- FindMarkers(pbmc, ident.1 = 1, ident.2 = 2)

# ============================================
# VISUALIZATION OF MARKER GENES
# ============================================

# Feature plots for top markers
top_genes <- top10_markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 2) %>% 
  pull(gene)

FeaturePlot(pbmc, features = top_genes[1:9], ncol = 3)

# Violin plots
VlnPlot(pbmc, features = top_genes[1:6], ncol = 3)

# Dot plot for marker genes
DotPlot(pbmc, features = unique(top10_markers$gene)) + 
  RotatedAxis()

# Heatmap of top markers
DoHeatmap(pbmc, features = top10_markers$gene) + 
  NoLegend()

# ============================================
# CELL TYPE ANNOTATION
# ============================================

# Define canonical PBMC markers
pbmc_markers <- list(
  "CD4+ T cells" = c("IL7R", "CD4", "TCF7", "LEF1"),
  "CD8+ T cells" = c("CD8A", "CD8B", "GZMK", "CCL5"),
  "Naive T cells" = c("CCR7", "SELL", "TCF7"),
  "Memory T cells" = c("IL7R", "S100A4"),
  "NK cells" = c("GNLY", "NKG7", "GZMB", "FGFBP2", "FCGR3A"),
  "B cells" = c("MS4A1", "CD79A", "CD79B", "IGHM"),
  "Monocytes (CD14+)" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Monocytes (CD16+)" = c("FCGR3A", "MS4A7", "CD16"),
  "Dendritic cells" = c("FCER1A", "CST3", "CLEC10A"),
  "Platelets" = c("PPBP", "PF4", "GNG11")
)

# Visualize canonical markers
for (cell_type in names(pbmc_markers)) {
  print(FeaturePlot(pbmc, features = pbmc_markers[[cell_type]], 
                    ncol = 2) + 
          plot_annotation(title = cell_type))
}

# Dot plot of canonical markers
DotPlot(pbmc, features = pbmc_markers) + 
  RotatedAxis() +
  ggtitle("Canonical PBMC Markers")

# ============================================
# AUTOMATED ANNOTATION
# ============================================

# Manual annotation based on marker expression
# Check average expression per cluster
avg_expression <- AverageExpression(pbmc, 
                                   features = unlist(pbmc_markers),
                                   return.seurat = FALSE)

# Assign cell type labels (example - adjust based on your data)
new_cluster_ids <- c(
  "0" = "CD4+ T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8+ T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "CD16+ Monocytes",
  "6" = "Dendritic cells",
  "7" = "Platelets"
)

# Rename clusters
pbmc$cell_type <- Idents(pbmc)
pbmc$cell_type <- plyr::mapvalues(
  x = pbmc$cell_type,
  from = names(new_cluster_ids),
  to = new_cluster_ids
)

Idents(pbmc) <- "cell_type"

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Annotated PBMC Cell Types")

# ============================================
# ALTERNATIVE: SingleR for automated annotation
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
pbmc$SingleR_labels <- predictions$labels

# Visualize SingleR predictions
DimPlot(pbmc, reduction = "umap", group.by = "SingleR_labels", 
        label = TRUE, pt.size = 0.5) +
  ggtitle("SingleR Automated Annotation")

# ============================================
# EXPORT RESULTS
# ============================================

# Save marker genes
write.csv(all_markers, "pbmc_all_markers.csv", row.names = FALSE)
write.csv(top10_markers, "pbmc_top10_markers.csv", row.names = FALSE)

# Save annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")

# Create summary table
cluster_summary <- data.frame(
  cluster = names(table(Idents(pbmc))),
  n_cells = as.numeric(table(Idents(pbmc))),
  cell_type = names(new_cluster_ids)
)

write.csv(cluster_summary, "cluster_summary.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
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

# Method 2: Find markers for specific cluster (e.g., cluster 0 vs all others)
cluster0_markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0_markers, n = 10)

# Method 3: Compare two specific clusters
cluster0_vs_1 <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1, min.pct = 0.25)

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
  "Naive CD4+ T" = c("IL7R", "CCR7", "CD4"),
  "Memory CD4+ T" = c("IL7R", "S100A4", "CD4"),
  "CD8+ T" = c("CD8A", "CD8B"),
  "NK" = c("GNLY", "NKG7", "NCAM1"),
  "B" = c("MS4A1", "CD79A"),
  "CD14+ Mono" = c("CD14", "LYZ"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
  "DC" = c("FCER1A", "CST3"),
  "Platelet" = c("PPBP")
)

# Check expression of canonical markers
for (cell_type in names(pbmc_markers)) {
  cat("\n", cell_type, ":\n")
  print(FeaturePlot(pbmc, features = pbmc_markers[[cell_type]][1:2]))
}

# Manual annotation based on marker expression
# Examine your top markers and compare with canonical markers
new_cluster_ids <- c(
  "0" = "CD14+ Mono",
  "1" = "CD4+ T",
  "2" = "CD8+ T",
  "3" = "B",
  "4" = "NK",
  "5" = "FCGR3A+ Mono",
  "6" = "DC",
  "7" = "Platelet"
)

# Rename clusters
pbmc$cell_type <- new_cluster_ids[as.character(Idents(pbmc))]
Idents(pbmc) <- "cell_type"

# Visualize annotated cell types
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Cell Type Annotation")

# ============================================
# AUTOMATED ANNOTATION (Optional)
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

# Compare manual vs automated annotation
table(pbmc$cell_type, pbmc$singler_labels)

# ============================================
# QUALITY CONTROL FOR MARKERS
# ============================================

# Check marker specificity
# Calculate percentage of cells expressing each marker in each cluster
marker_expression <- function(seurat_obj, markers) {
  results <- data.frame()
  for (marker in markers) {
    pct <- sapply(levels(Idents(seurat_obj)), function(cluster) {
      cells <- WhichCells(seurat_obj, idents = cluster)
      sum(GetAssayData(seurat_obj, slot = "data")[marker, cells] > 0) / length(cells) * 100
    })
    results <- rbind(results, data.frame(Gene = marker, Cluster = names(pct), Percentage = pct))
  }
  return(results)
}

# Example usage
key_markers <- c("CD3D", "CD14", "MS4A1", "GNLY")
expr_summary <- marker_expression(pbmc, key_markers)
print(expr_summary)

# ============================================
# EXPORT RESULTS
# ============================================

# Save annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")

# Export cell type assignments
cell_annotations <- data.frame(
  Cell = colnames(pbmc),
  Cluster = Idents(pbmc),
  CellType = pbmc$cell_type,
  UMAP1 = Embeddings(pbmc, "umap")[,1],
  UMAP2 = Embeddings(pbmc, "umap")[,2]
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)

# Create summary statistics
cluster_summary <- table(pbmc$cell_type)
print(cluster_summary)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
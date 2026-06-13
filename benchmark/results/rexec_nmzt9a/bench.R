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

# ============================================================================
# 1. DATA LOADING AND QC
# ============================================================================

# Load 10X data (adjust path to your data)
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "PBMC_Analysis", 
                           min.cells = 3, 
                           min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

# ============================================================================
# 2. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# ============================================================================
# 3. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 4. SCALING AND PCA
# ============================================================================

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 5. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 6. FIND CLUSTER MARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top_markers)

# Visualize top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top5$gene) + NoLegend()

# ============================================================================
# 7. IMMUNE CELL TYPE ANNOTATION
# ============================================================================

# Define canonical immune cell markers
immune_markers <- list(
  "CD4_T" = c("CD3D", "CD3E", "CD4", "IL7R"),
  "CD8_T" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "NK" = c("GNLY", "NKG7", "NCAM1", "KLRD1"),
  "B_cells" = c("MS4A1", "CD79A", "CD79B"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "DC" = c("FCER1A", "CST3", "CD1C"),
  "Platelets" = c("PPBP", "PF4")
)

# Visualize marker expression
for (cell_type in names(immune_markers)) {
  print(FeaturePlot(pbmc, features = immune_markers[[cell_type]], 
                    ncol = 2))
}

# Dot plot of canonical markers
DotPlot(pbmc, features = unlist(immune_markers)) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plots for key markers
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY"), 
        ncol = 2)

# ============================================================================
# 8. MANUAL ANNOTATION (adjust based on your marker expression)
# ============================================================================

# Example annotation - adjust cluster numbers based on your results
new.cluster.ids <- c(
  "CD4 T cells",      # 0
  "CD14+ Monocytes",  # 1
  "B cells",          # 2
  "CD8 T cells",      # 3
  "NK cells",         # 4
  "FCGR3A+ Monocytes",# 5
  "Dendritic cells",  # 6
  "Platelets"         # 7
)

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# ============================================================================
# 9. DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# ============================================================================

# Example: Compare CD4 vs CD8 T cells
cd4_vs_cd8 <- FindMarkers(pbmc, 
                          ident.1 = "CD4 T cells", 
                          ident.2 = "CD8 T cells")

head(cd4_vs_cd8, n = 20)

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_analyzed.rds")

# Export marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "top_markers_per_cluster.csv", row.names = FALSE)

# Export cell type annotations
cell_annotations <- data.frame(
  cell_barcode = colnames(pbmc),
  cell_type = Idents(pbmc),
  cluster = pbmc@meta.data$seurat_clusters
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)

# ============================================================================
# 11. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Feature plots for key markers
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY", 
                               "CD8A", "CD4", "FCGR3A", "FCER1A"))

# Split by cell type
DimPlot(pbmc, reduction = "umap", split.by = "ident", ncol = 4)

# Cell type proportions
cell_counts <- table(Idents(pbmc))
cell_props <- prop.table(cell_counts) * 100

barplot(cell_props, 
        las = 2, 
        main = "Cell Type Proportions",
        ylab = "Percentage",
        col = rainbow(length(cell_props)))

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
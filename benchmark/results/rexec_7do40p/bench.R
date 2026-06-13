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
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC_Analysis",
  min.cells = 3,      # Keep genes expressed in >= 3 cells
  min.features = 200  # Keep cells with >= 200 genes
)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)

# Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = 
  nFeature_RNA > 200 & 
  nFeature_RNA < 5000 & 
  percent.mt < 15
)

# ============================================================================
# 2. NORMALIZATION AND FEATURE SELECTION
# ============================================================================

# Normalize data (log-normalization)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 3. SCALING AND PCA
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
# 4. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 5. FIND CLUSTER MARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(
  pbmc, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top_markers)

# Save marker genes
write.csv(pbmc.markers, "pbmc_cluster_markers.csv", row.names = FALSE)

# ============================================================================
# 6. VISUALIZE KEY IMMUNE CELL MARKERS
# ============================================================================

# Define canonical immune cell markers
immune_markers <- list(
  "T cells (CD3)" = c("CD3D", "CD3E", "CD3G"),
  "CD4+ T cells" = c("CD4", "IL7R"),
  "CD8+ T cells" = c("CD8A", "CD8B"),
  "NK cells" = c("GNLY", "NKG7", "NCAM1"),
  "B cells" = c("MS4A1", "CD79A", "CD79B"),
  "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Dendritic cells" = c("FCER1A", "CST3"),
  "Platelets" = c("PPBP")
)

# Feature plots for key markers
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY", 
                               "CD8A", "FCGR3A", "LYZ", "PPBP"), 
            ncol = 4)

# Dot plot of markers
DotPlot(pbmc, features = unlist(immune_markers)) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plots for specific markers
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY"), 
        ncol = 2, pt.size = 0)

# Heatmap of top markers
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + 
  NoLegend()

# ============================================================================
# 7. CELL TYPE ANNOTATION
# ============================================================================

# Manual annotation based on marker expression
new.cluster.ids <- c(
  "CD14+ Monocytes",  # 0
  "CD4+ T cells",     # 1
  "CD8+ T cells",     # 2
  "B cells",          # 3
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

# Save annotated object
pbmc$cell_type <- Idents(pbmc)

# ============================================================================
# 8. CELL TYPE PROPORTIONS
# ============================================================================

# Calculate cell type proportions
cell_counts <- table(Idents(pbmc))
cell_proportions <- prop.table(cell_counts) * 100

# Visualize proportions
barplot(cell_proportions, 
        las = 2, 
        main = "Cell Type Proportions",
        ylab = "Percentage (%)",
        col = rainbow(length(cell_proportions)))

# ============================================================================
# 9. DIFFERENTIAL EXPRESSION BETWEEN CELL TYPES
# ============================================================================

# Example: Compare CD4+ vs CD8+ T cells
cd4_vs_cd8 <- FindMarkers(
  pbmc, 
  ident.1 = "CD4+ T cells", 
  ident.2 = "CD8+ T cells",
  min.pct = 0.25
)

head(cd4_vs_cd8, n = 20)

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_analyzed.rds")

# Export UMAP coordinates with cell types
umap_coords <- as.data.frame(Embeddings(pbmc, reduction = "umap"))
umap_coords$cell_type <- Idents(pbmc)
write.csv(umap_coords, "umap_coordinates_with_celltypes.csv")

# Export expression matrix
write.csv(as.matrix(GetAssayData(pbmc, slot = "data")), 
          "normalized_expression_matrix.csv")

# ============================================================================
# 11. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Split UMAP by cell type
DimPlot(pbmc, reduction = "umap", split.by = "cell_type", ncol = 4)

# Ridge plot for marker genes
RidgePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "GNLY"), ncol = 2)

# Generate summary statistics
cat("\n=== Analysis Summary ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Number of cell types identified:", length(unique(Idents(pbmc))), "\n")
cat("\nCell type distribution:\n")
print(cell_counts)
cat("\nCell type percentages:\n")
print(round(cell_proportions, 2))

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
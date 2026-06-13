setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_analysis", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC visualization
qc_vln <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("qc_violin_plot.png", plot = qc_vln, width = 12, height = 4)

qc_scatter1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc_scatter2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc_combined <- qc_scatter1 + qc_scatter2
ggsave("qc_scatter_plots.png", plot = qc_combined, width = 12, height = 5)

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify top 10 variable features
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features
var_plot1 <- VariableFeaturePlot(pbmc)
var_plot2 <- LabelPoints(plot = var_plot1, points = top10, repel = TRUE)
ggsave("variable_features.png", plot = var_plot2, width = 10, height = 6)

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA
pca_plot <- DimPlot(pbmc, reduction = "pca")
ggsave("pca_plot.png", plot = pca_plot, width = 8, height = 6)

elbow_plot <- ElbowPlot(pbmc)
ggsave("elbow_plot.png", plot = elbow_plot, width = 8, height = 6)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
set.seed(42)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("umap_clusters.png", plot = umap_plot, width = 10, height = 8)

# Find all markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(pbmc.markers, "all_cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "top10_markers_per_cluster.csv", row.names = FALSE)

# Visualize top markers
top5 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

heatmap_plot <- DoHeatmap(pbmc, features = top5$gene) + NoLegend()
ggsave("marker_heatmap.png", plot = heatmap_plot, width = 12, height = 10)

# Feature plots for selected markers
if(nrow(top10) > 0) {
  feature_genes <- head(top10, 9)
  feature_plot <- FeaturePlot(pbmc, features = feature_genes)
  ggsave("feature_plot_top_genes.png", plot = feature_plot, width = 15, height = 12)
}

# Violin plots for selected markers
if(nrow(top10) > 0) {
  violin_genes <- head(top10, 9)
  violin_plot <- VlnPlot(pbmc, features = violin_genes)
  ggsave("violin_plot_top_genes.png", plot = violin_plot, width = 15, height = 12)
}

# Save Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

print("Analysis complete. All outputs saved.")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
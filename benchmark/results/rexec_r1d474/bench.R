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
library(patchwork)

# Set working directory and output folder
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read raw count matrix (genes x cells)
# Replace 'counts_matrix.csv' with your actual file path
counts <- read.csv("counts_matrix.csv", row.names = 1, check.names = FALSE)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
pdf(file.path(output_dir, "QC_violin_plots.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data using SCTransform
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Run PCA
pbmc <- RunPCA(pbmc, verbose = FALSE)

# Determine dimensionality
pdf(file.path(output_dir, "ElbowPlot.pdf"))
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)

# Save Seurat object
saveRDS(pbmc, file = file.path(output_dir, "pbmc_seurat_object.rds"))
SaveH5Seurat(pbmc, filename = file.path(output_dir, "pbmc_seurat_object.h5Seurat"))

# Plot UMAP colored by cluster
umap_cluster_plot <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + ggtitle("UMAP - Clusters")
ggsave(filename = file.path(output_dir, "UMAP_clusters.png"), plot = umap_cluster_plot, width = 8, height = 6, dpi = 300)

# Plot UMAP colored by percent.mt
umap_mt_plot <- DimPlot(pbmc, reduction = "umap", group.by = "percent.mt", pt.size = 0.5) + ggtitle("UMAP - %MT")
ggsave(filename = file.path(output_dir, "UMAP_percentMT.png"), plot = umap_mt_plot, width = 8, height = 6, dpi = 300)

# Identify cluster markers
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = file.path(output_dir, "cluster_markers.csv"), row.names = FALSE)

# Save top markers per cluster
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top_markers, file = file.path(output_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)

# Save feature plots for selected markers
selected_genes <- c("CD3D", "CD14", "MS4A1", "GNLY", "NKG7")
for (gene in selected_genes) {
  p <- FeaturePlot(pbmc, features = gene, reduction = "umap") + ggtitle(paste("UMAP -", gene))
  ggsave(filename = file.path(output_dir, paste0("FeaturePlot_", gene, ".png")), plot = p, width = 6, height = 5, dpi = 300)
}

# Save dimensionality reduction plots
pdf(file.path(output_dir, "PCA_plot.pdf"))
DimPlot(pbmc, reduction = "pca", dims = c(1,2), pt.size = 0.5) + ggtitle("PCA - PC1 vs PC2")
dev.off()

# Save cluster heatmap of top markers
top10_genes <- top_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% pull(gene) %>% unique()
DoHeatmap(pbmc, features = top10_genes, group.by = "seurat_clusters", size = 3) + ggtitle("Heatmap of Top Markers") %>% 
  ggsave(filename = file.path(output_dir, "Heatmap_top_markers.png"), width = 10, height = 8, dpi = 300)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
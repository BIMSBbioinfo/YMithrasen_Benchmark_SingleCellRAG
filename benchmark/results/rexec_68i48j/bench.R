setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Set random seed for reproducibility
set.seed(1234)

# Define file paths
raw_counts_path <- "data/raw_counts"          # Path to 10x Genomics raw counts folder
output_dir <- "analysis_results"
dir.create(output_dir, showWarnings = FALSE)

# 1. Read raw count data
counts <- Read10X(data.dir = raw_counts_path)

# 2. Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# 3. Add mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4. Visualize QC metrics
pdf(file = file.path(output_dir, "QC_violin_plots.pdf"))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 5. Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 6. Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 7. Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 8. Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 9. Perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 10. Determine dimensionality
pdf(file = file.path(output_dir, "ElbowPlot.pdf"))
ElbowPlot(pbmc, ndims = 50)
dev.off()

# 11. Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 12. Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# 13. Visualize UMAP
pdf(file = file.path(output_dir, "UMAP_clusters.pdf"))
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of PBMCs")
dev.off()

# 14. Identify cluster markers
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = file.path(output_dir, "cluster_markers.csv"), row.names = FALSE)

# 15. Save Seurat object
saveRDS(pbmc, file = file.path(output_dir, "pbmc_seurat_object.rds"))
SaveH5Seurat(pbmc, filename = file.path(output_dir, "pbmc_seurat_object.h5Seurat"))

# 16. Save dimensionality reduction plots
pdf(file = file.path(output_dir, "PCA_plot.pdf"))
DimPlot(pbmc, reduction = "pca", group.by = "seurat_clusters") + ggtitle("PCA of PBMCs")
dev.off()

# 17. Save feature plots for top markers of each cluster
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
for (clust in unique(top_markers$cluster)) {
  markers <- top_markers %>% filter(cluster == clust) %>% pull(gene)
  pdf(file = file.path(output_dir, paste0("FeaturePlot_cluster_", clust, ".pdf")))
  FeaturePlot(pbmc, features = markers, cols = c("lightgrey", "blue")) + ggtitle(paste("Cluster", clust, "Markers"))
  dev.off()
}
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
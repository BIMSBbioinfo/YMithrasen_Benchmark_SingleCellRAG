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

set.seed(42)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = ".")

# Initialize the Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_analysis", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC visualization
pdf("qc_plots.pdf", width = 12, height = 4)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("feature_scatter_plots.pdf", width = 10, height = 5)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization using SCTransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

pdf("pca_elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30)

pdf("umap_clusters.pdf", width = 8, height = 6)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# Find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save top markers
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(pbmc.markers, "all_markers.csv", row.names = FALSE)
write.csv(top10_markers, "top10_markers.csv", row.names = FALSE)

# Visualize top markers
pdf("top_markers_heatmap.pdf", width = 12, height = 10)
DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()
dev.off()

# Save Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
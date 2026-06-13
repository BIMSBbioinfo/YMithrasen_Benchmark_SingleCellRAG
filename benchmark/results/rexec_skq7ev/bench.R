setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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
sc.data <- Read10X(data.dir = ".")

# Initialize the Seurat object
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc", min.cells = 3, min.features = 200)

# Add mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# QC visualization
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Visualize top 10 variable features
top10 <- head(VariableFeatures(sc), 10)
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1
plot2

# Scale data
sc <- ScaleData(sc, features = VariableFeatures(sc))

# Run PCA
sc <- RunPCA(sc, features = VariableFeatures(sc))

# Visualize PCA
DimPlot(sc, reduction = "pca")
ElbowPlot(sc)

# Clustering
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

# Run UMAP
sc <- RunUMAP(sc, dims = 1:10)

# Visualize UMAP
DimPlot(sc, reduction = "umap", label = TRUE)

# Find all markers
sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Visualize canonical markers
FeaturePlot(sc, features = c("CD3D", "CD3E", "CD14", "LYZ", "CD79A", "MS4A1", "NKG7", "GNLY"))

VlnPlot(sc, features = c("CD3D", "CD14", "CD79A", "NKG7"))

# Heatmap of top markers
sc.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_markers
DoHeatmap(sc, features = top10_markers$gene) + NoLegend()

# Assign cell type identities
new.cluster.ids <- c("CD4 T", "CD14 Mono", "B", "CD8 T", "NK", "FCGR3A Mono", "DC", "Platelet")
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)

# Final UMAP with cell types
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
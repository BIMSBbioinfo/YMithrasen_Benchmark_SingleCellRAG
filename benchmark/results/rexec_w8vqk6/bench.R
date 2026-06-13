setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes and genes
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes$V2  # Use HGNC gene symbols
colnames(matrix) <- barcodes$V1  # Use cell barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Loaded data with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 8. UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of PBMC Clusters")

# ============================================================================
# 9. FIND MARKER GENES FOR ALL CLUSTERS
# ============================================================================

# Find markers for every cluster compared to all remaining cells
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# View top markers
cat("\nTop 5 markers per cluster:\n")
top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
print(top5_markers)

# Save all markers to file
write.csv(pbmc.markers, "pbmc_all_cluster_markers.csv", row.names = FALSE)

# ============================================================================
# 10. HEATMAP OF TOP 10 MARKERS PER CLUSTER
# ============================================================================

# Extract top 10 markers per cluster
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# Create heatmap
DoHeatmap(pbmc, features = top10_markers$gene, size = 3) + 
  ggtitle("Top 10 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 6))

# Alternative: More customized heatmap
DoHeatmap(pbmc, 
          features = top10_markers$gene, 
          group.by = "seurat_clusters",
          size = 3,
          angle = 0,
          hjust = 0.5) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  ggtitle("Top 10 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 6))

# ============================================================================
# 11. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Feature plots for top markers
top_genes <- top10_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 1, order_by = avg_log2FC) %>% 
  pull(gene)

FeaturePlot(pbmc, features = top_genes[1:min(9, length(top_genes))], ncol = 3)

# Violin plots for top markers
VlnPlot(pbmc, features = top_genes[1:min(6, length(top_genes))], ncol = 3)

# Dot plot
DotPlot(pbmc, features = unique(top10_markers$gene)) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 6))

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# Save UMAP plot
ggsave("pbmc_umap_clusters.pdf", width = 8, height = 6)

cat("\nAnalysis complete! Results saved.\n")
cat("Seurat object saved as: pbmc_seurat_object.rds\n")
cat("Marker genes saved as: pbmc_all_cluster_markers.csv\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
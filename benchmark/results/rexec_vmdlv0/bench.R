setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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

print(paste("Loaded", ncol(pbmc), "cells and", nrow(pbmc), "genes"))

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

print(paste("After QC:", ncol(pbmc), "cells retained"))

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
print("Top 10 variable genes:")
print(top10)

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

# Check cluster assignments
print(table(Idents(pbmc)))

# ============================================================================
# 8. UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC UMAP Clustering")

# ============================================================================
# 9. FIND MARKER GENES FOR ALL CLUSTERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# View top markers
print(pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC))

# Save all markers to file
write.csv(pbmc.markers, "pbmc_all_cluster_markers.csv", row.names = FALSE)

# ============================================================================
# 10. SELECT TOP 10 MARKERS PER CLUSTER FOR HEATMAP
# ============================================================================

# Get top 10 markers per cluster
top10_markers <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

print(paste("Selected", nrow(top10_markers), "marker genes for heatmap"))

# ============================================================================
# 11. CREATE HEATMAP
# ============================================================================

# Create heatmap with DoHeatmap
heatmap <- DoHeatmap(pbmc, features = top10_markers$gene) + 
  ggtitle("Top 10 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 6))

print(heatmap)

# Save heatmap
ggsave("pbmc_top10_markers_heatmap.pdf", plot = heatmap, width = 12, height = 16)

# ============================================================================
# 12. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Feature plots for known cell type markers
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Violin plots for top markers
VlnPlot(pbmc, features = head(top10_markers$gene, 9), ncol = 3)

# Dot plot for top markers
top5_markers <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

DotPlot(pbmc, features = unique(top5_markers$gene)) + 
  RotatedAxis() +
  ggtitle("Top 5 Markers per Cluster")

# ============================================================================
# 13. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# Save UMAP coordinates
umap_coords <- as.data.frame(Embeddings(pbmc, reduction = "umap"))
umap_coords$cluster <- Idents(pbmc)
write.csv(umap_coords, "pbmc_umap_coordinates.csv")

print("Analysis complete! All results saved.")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
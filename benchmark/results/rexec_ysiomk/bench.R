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
library(ggplot2)

# ============================================
# 1. LOAD DATA
# ============================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes (cell IDs)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene IDs and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Loaded data with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================
# 2. QUALITY CONTROL & PREPROCESSING
# ============================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After filtering:", ncol(pbmc), "cells remaining\n")

# ============================================
# 3. NORMALIZATION
# ============================================

# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================
# 4. FEATURE SELECTION
# ============================================

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================
# 5. SCALING
# ============================================

# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================
# 6. PCA
# ============================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================
# 7. CLUSTERING
# ============================================

# Find neighbors using first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================
# 8. UMAP EMBEDDING
# ============================================

# Run UMAP using first 10 PCs
pbmc <- RunUMAP(pbmc, dims = 1:10)

# ============================================
# 9. VISUALIZATION
# ============================================

# Plot UMAP with clusters
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("PBMC Clusters on UMAP") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# Save the plot
ggsave("pbmc_umap_clusters.png", plot = p1, width = 10, height = 8, dpi = 300)

# Additional visualization: color by number of genes
p2 <- FeaturePlot(pbmc, reduction = "umap", features = "nFeature_RNA")
print(p2)

# ============================================
# 10. CLUSTER STATISTICS
# ============================================

# Print cluster sizes
cluster_sizes <- table(Idents(pbmc))
cat("\nCluster sizes:\n")
print(cluster_sizes)

# ============================================
# 11. SAVE RESULTS
# ============================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analyzed.rds")

cat("\nAnalysis complete! Results saved.\n")

# Optional: Find marker genes for each cluster
# This helps identify cell types
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

cat("\nTop 5 marker genes per cluster:\n")
print(top5_markers)

# Visualize top markers
top_genes <- unique(top5_markers$gene)[1:9]  # Take first 9 unique genes
FeaturePlot(pbmc, features = top_genes, ncol = 3)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
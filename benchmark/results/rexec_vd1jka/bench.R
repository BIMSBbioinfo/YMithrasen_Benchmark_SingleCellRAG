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
library(Matrix)
library(ggplot2)

# ============================================================================
# 1. DATA LOADING
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
                           min.cells = 3,      # Keep genes expressed in >= 3 cells
                           min.features = 200) # Keep cells with >= 200 genes

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
print(paste("Top 10 variable genes:", paste(top10, collapse = ", ")))

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
# 6. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality of the dataset
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check number of cells per cluster
table(Idents(pbmc))

# ============================================================================
# 8. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP)
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of PBMC 2,700 cells")

# Save the plot
ggsave("UMAP_clusters.pdf", width = 8, height = 6)

# ============================================================================
# 9. FIND ALL MARKER GENES
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE,           # Only positive markers
                               min.pct = 0.25,            # Min % of cells expressing
                               logfc.threshold = 0.25)    # Log2FC threshold

# View top markers
print(head(pbmc.markers, n = 20))

# Get top 10 markers per cluster
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top10_markers)

# Save marker genes to file
write.csv(pbmc.markers, "all_cluster_markers.csv", row.names = FALSE)
write.csv(top10_markers, "top10_markers_per_cluster.csv", row.names = FALSE)

# ============================================================================
# 10. VISUALIZATION OF MARKER GENES
# ============================================================================

# Heatmap of top markers
DoHeatmap(pbmc, features = top10_markers$gene) + 
  theme(axis.text.y = element_text(size = 6))
ggsave("marker_heatmap.pdf", width = 12, height = 10)

# Violin plots for top marker of each cluster
top_marker_per_cluster <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC)

VlnPlot(pbmc, features = top_marker_per_cluster$gene[1:min(9, nrow(top_marker_per_cluster))])

# Feature plots for selected markers
FeaturePlot(pbmc, features = top_marker_per_cluster$gene[1:min(9, nrow(top_marker_per_cluster))])

# ============================================================================
# 11. SAVE SEURAT OBJECT
# ============================================================================

saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells after QC:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Total marker genes found:", nrow(pbmc.markers), "\n")
cat("\nCells per cluster:\n")
print(table(Idents(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
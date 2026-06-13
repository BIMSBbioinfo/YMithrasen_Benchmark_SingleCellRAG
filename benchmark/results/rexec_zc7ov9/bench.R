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

# ============================================================================
# 2. CREATE SEURAT OBJECT
# ============================================================================

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3,      # Keep genes expressed in >= 3 cells
                           min.features = 200) # Keep cells with >= 200 genes

# ============================================================================
# 3. QUALITY CONTROL
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
# 4. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# ============================================================================
# 5. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
print(paste("Top 10 variable genes:", paste(top10, collapse = ", ")))

# ============================================================================
# 6. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 7. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 8. CLUSTERING
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ============================================================================
# 9. UMAP VISUALIZATION
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 10. FIND MARKER GENES FOR ALL CLUSTERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE,      # Only positive markers
                               min.pct = 0.25,       # Min % of cells expressing
                               logfc.threshold = 0.25)

# View top markers
print(pbmc.markers %>% 
      group_by(cluster) %>% 
      top_n(n = 5, wt = avg_log2FC))

# ============================================================================
# 11. SELECT TOP 10 MARKERS PER CLUSTER
# ============================================================================

# Get top 10 markers per cluster
top10_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# ============================================================================
# 12. HEATMAP VISUALIZATION
# ============================================================================

# Create heatmap of top 10 markers per cluster
DoHeatmap(pbmc, features = top10_markers$gene) + 
  NoLegend() +
  theme(axis.text.y = element_text(size = 6))

# Alternative: with legend and better formatting
DoHeatmap(pbmc, features = top10_markers$gene, size = 3) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.y = element_text(size = 8))

# ============================================================================
# 13. SAVE RESULTS
# ============================================================================

# Save marker genes to CSV
write.csv(pbmc.markers, "pbmc_all_markers.csv", row.names = FALSE)
write.csv(top10_markers, "pbmc_top10_markers_per_cluster.csv", row.names = FALSE)

# Save plots
ggsave("umap_clusters.png", DimPlot(pbmc, reduction = "umap", label = TRUE), 
       width = 10, height = 8)

pdf("heatmap_top10_markers.pdf", width = 12, height = 10)
DoHeatmap(pbmc, features = top10_markers$gene, size = 3) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

# Save Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# ============================================================================
# 14. SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells after QC:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(pbmc)))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
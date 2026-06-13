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
library(dplyr)

# ============================================================================
# 1. Load the 10X-style data files
# ============================================================================

# Read the sparse matrix (Matrix Market format)
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene IDs and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# ============================================================================
# 2. Create Seurat object
# ============================================================================

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Keep genes expressed in at least 3 cells
  min.features = 200  # Keep cells with at least 200 genes
)

# ============================================================================
# 3. Quality Control (QC) filtering
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("Cells after QC filtering:", ncol(pbmc), "\n")

# ============================================================================
# 4. Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 5. Identify highly variable features
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================================================
# 6. Scaling the data
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 7. Perform linear dimensional reduction (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality of the dataset
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 8. Clustering
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 9. Run UMAP for visualization
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC 2700 - UMAP Clustering")

# ============================================================================
# 10. Find marker genes for every cluster
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(
  pbmc, 
  only.pos = TRUE,           # Only positive markers
  min.pct = 0.25,            # Expressed in at least 25% of cells
  logfc.threshold = 0.25     # Log fold-change threshold
)

# View top markers
cat("\nTop markers per cluster:\n")
print(pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))

# ============================================================================
# 11. Visualize top 5 markers per cluster as a heatmap
# ============================================================================

# Get top 5 markers per cluster
top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Create heatmap
DoHeatmap(
  pbmc, 
  features = top5_markers$gene,
  size = 3,
  angle = 90
) + 
  ggtitle("Top 5 Marker Genes per Cluster") +
  theme(axis.text.y = element_text(size = 6))

# ============================================================================
# 12. Save results
# ============================================================================

# Save marker genes to CSV
write.csv(pbmc.markers, "pbmc_all_markers.csv", row.names = FALSE)
write.csv(top5_markers, "pbmc_top5_markers_per_cluster.csv", row.names = FALSE)

# Save Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

cat("\nAnalysis complete! Results saved.\n")
cat("- Marker genes: pbmc_all_markers.csv\n")
cat("- Top 5 markers: pbmc_top5_markers_per_cluster.csv\n")
cat("- Seurat object: pbmc_seurat_object.rds\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
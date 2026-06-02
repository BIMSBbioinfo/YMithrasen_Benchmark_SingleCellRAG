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

# Set working directory (adjust to your file location)
# setwd(".")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (feature identifiers)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_name")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_name
colnames(matrix) <- barcodes$barcode

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 2. QUALITY CONTROL (Optional but recommended)
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics (adjust thresholds as needed)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

# Normalize the data using log-normalization
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

cat("Data normalized\n")

# ============================================================================
# 4. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================

# Find the top 2000 highly variable genes
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# 5. SCALE THE DATA
# ============================================================================

# Scale all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaled\n")

# ============================================================================
# 6. PERFORM PCA
# ============================================================================

# Run PCA on the scaled data
pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc), 
               npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# Visualize PCA
DimPlot(pbmc, reduction = "pca")

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
# Using first 20 PCs (adjust based on elbow plot)
pbmc <- FindNeighbors(pbmc, dims = 1:20)

# Find clusters with resolution 0.5 (adjust for granularity)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 8. RUN UMAP
# ============================================================================

# Compute UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

cat("UMAP computed\n")

# ============================================================================
# 9. VISUALIZE CLUSTERS
# ============================================================================

# UMAP plot colored by cluster
p1 <- DimPlot(pbmc, 
              reduction = "umap", 
              label = TRUE, 
              label.size = 6,
              pt.size = 1) +
  ggtitle("PBMC Clusters (UMAP)")
print(p1)

# Save high-resolution plot
ggsave("pbmc_umap_clusters.png", plot = p1, width = 10, height = 8, dpi = 300)

# UMAP plot with cluster labels
p2 <- DimPlot(pbmc, 
              reduction = "umap", 
              label = TRUE, 
              label.size = 5,
              repel = TRUE,
              pt.size = 1) +
  NoLegend() +
  ggtitle("PBMC Clusters with Labels")
print(p2)

# ============================================================================
# 10. FIND CLUSTER MARKERS (Optional)
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Get top 5 markers per cluster
top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top5_markers)

# Visualize top markers
top_genes <- top5_markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 2) %>% 
  pull(gene) %>% 
  unique()

FeaturePlot(pbmc, features = head(top_genes, 9), ncol = 3)

# Heatmap of top markers
DoHeatmap(pbmc, features = top5_markers$gene) + 
  NoLegend() +
  theme(text = element_text(size = 8))

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Save cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

# Save marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")
cat("Seurat object saved as: pbmc_analysis.rds\n")
cat("Cluster assignments saved as: cluster_assignments.csv\n")
cat("Marker genes saved as: cluster_markers.csv\n")

# ============================================================================
# 12. SUMMARY STATISTICS
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
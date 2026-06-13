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

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
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

cat("Loaded", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC:", ncol(pbmc), "cells remaining\n")

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

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 8. UMAP VISUALIZATION
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot UMAP with cluster labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Clusters (Unannotated)")

# ============================================================================
# 9. FIND CLUSTER MARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Display top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# ============================================================================
# 10. VISUALIZE CANONICAL MARKERS
# ============================================================================

# Define canonical markers for PBMC cell types
canonical_markers <- c(
  "IL7R", "CCR7",    # CD4 T cells
  "CD14", "LYZ",     # CD14 Monocytes
  "MS4A1",           # B cells
  "CD8A",            # CD8 T cells
  "GNLY", "NKG7"     # NK cells
)

# Feature plots for canonical markers
FeaturePlot(pbmc, features = canonical_markers, ncol = 3)

# Dot plot for canonical markers
DotPlot(pbmc, features = canonical_markers) + 
  RotatedAxis() +
  ggtitle("Canonical Marker Expression by Cluster")

# Violin plots for key markers
VlnPlot(pbmc, features = canonical_markers, ncol = 3)

# ============================================================================
# 11. CELL TYPE ANNOTATION
# ============================================================================

# Create a named vector for cluster annotation
# Adjust cluster numbers based on your marker expression patterns
# This is a template - you'll need to adjust based on your actual results

# First, examine average expression of markers per cluster
avg_expression <- AverageExpression(pbmc, features = canonical_markers, group.by = "seurat_clusters")
print(avg_expression$RNA)

# Manual annotation based on marker expression
# IMPORTANT: Adjust these cluster assignments based on your marker visualization
new.cluster.ids <- c(
  "CD4 T cells",      # Cluster 0 (example - adjust based on IL7R, CCR7)
  "CD14 Monocytes",   # Cluster 1 (example - adjust based on CD14, LYZ)
  "CD4 T cells",      # Cluster 2
  "B cells",          # Cluster 3 (example - adjust based on MS4A1)
  "CD8 T cells",      # Cluster 4 (example - adjust based on CD8A)
  "NK cells",         # Cluster 5 (example - adjust based on GNLY, NKG7)
  "CD14 Monocytes",   # Cluster 6
  "CD4 T cells"       # Cluster 7
)

# Apply annotations
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Store cell type in metadata
pbmc$cell_type <- Idents(pbmc)

# ============================================================================
# 12. FINAL ANNOTATED UMAP
# ============================================================================

# Create annotated UMAP
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4) + 
  ggtitle("Annotated PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# Save high-resolution plot
ggsave("pbmc_annotated_umap.png", plot = p1, width = 10, height = 8, dpi = 300)

# ============================================================================
# 13. SUMMARY STATISTICS
# ============================================================================

# Cell type counts
cell_type_counts <- table(pbmc$cell_type)
print("Cell type distribution:")
print(cell_type_counts)
print(prop.table(cell_type_counts))

# ============================================================================
# 14. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")

# Export cell type annotations
annotations <- data.frame(
  barcode = colnames(pbmc),
  cluster = pbmc$seurat_clusters,
  cell_type = pbmc$cell_type
)
write.csv(annotations, "cell_type_annotations.csv", row.names = FALSE)

# Export marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)

cat("\nAnalysis complete! Files saved:\n")
cat("- pbmc_annotated_umap.png\n")
cat("- pbmc_annotated.rds\n")
cat("- cell_type_annotations.csv\n")
cat("- cluster_markers.csv\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
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
# STEP 1: Load the Data
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
# STEP 2: Quality Control
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 3: Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# STEP 4: Feature Selection
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================================================
# STEP 5: Scaling
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 6: PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# STEP 7: Clustering
# ============================================================================

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ============================================================================
# STEP 8: UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Clusters (Unnannotated)")

# ============================================================================
# STEP 9: Find Cluster Markers
# ============================================================================

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top_markers)

# ============================================================================
# STEP 10: Visualize Known Cell Type Markers
# ============================================================================

# Define marker genes for each cell type
marker_genes <- c(
  "IL7R", "CCR7",      # CD4 T cells
  "CD14", "LYZ",       # CD14 Monocytes
  "MS4A1",             # B cells
  "CD8A",              # CD8 T cells
  "GNLY", "NKG7"       # NK cells
)

# Feature plots for marker genes
FeaturePlot(pbmc, features = marker_genes, ncol = 3)

# Dot plot for marker expression across clusters
DotPlot(pbmc, features = marker_genes) + 
  RotatedAxis() + 
  ggtitle("Marker Gene Expression Across Clusters")

# Violin plots
VlnPlot(pbmc, features = marker_genes, ncol = 3)

# ============================================================================
# STEP 11: Cell Type Annotation
# ============================================================================

# Examine average expression of markers per cluster
avg_exp <- AverageExpression(pbmc, features = marker_genes, group.by = "seurat_clusters")
print(avg_exp$RNA)

# Manual annotation based on marker expression
# Adjust cluster numbers based on your actual results
new.cluster.ids <- c(
  "CD4 T cells",       # Cluster 0 (IL7R+, CCR7+)
  "CD14 Monocytes",    # Cluster 1 (CD14+, LYZ+)
  "CD4 T cells",       # Cluster 2 (IL7R+)
  "B cells",           # Cluster 3 (MS4A1+)
  "CD8 T cells",       # Cluster 4 (CD8A+)
  "NK cells",          # Cluster 5 (GNLY+, NKG7+)
  "CD14 Monocytes",    # Cluster 6 (CD14+, LYZ+)
  "CD4 T cells"        # Cluster 7 (IL7R+)
)

# Note: You'll need to adjust the above based on your actual cluster-marker associations
# Use the DotPlot and average expression to determine which cluster corresponds to which cell type

# Rename clusters
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Add cell type to metadata
pbmc$cell_type <- Idents(pbmc)

# ============================================================================
# STEP 12: Final Annotated UMAP
# ============================================================================

# Create annotated UMAP
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 4) + 
  ggtitle("Annotated PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# Save the plot
ggsave("annotated_pbmc_umap.png", plot = p1, width = 10, height = 8, dpi = 300)

# ============================================================================
# STEP 13: Summary Statistics
# ============================================================================

# Cell type counts
cell_type_counts <- table(pbmc$cell_type)
print("Cell Type Distribution:")
print(cell_type_counts)

# Proportions
print("Cell Type Proportions:")
print(prop.table(cell_type_counts))

# ============================================================================
# STEP 14: Save Results
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")

# Export cell type annotations
annotations <- data.frame(
  barcode = colnames(pbmc),
  cell_type = pbmc$cell_type,
  cluster = pbmc$seurat_clusters
)
write.csv(annotations, "cell_type_annotations.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Comprehensive Seurat Pipeline for PBMC Single-Cell Analysis

Here's a complete, executable R script to analyze your 2,700 PBMC single-cell dataset:

```r
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data
# ============================================================================

# Read the Matrix Market format sparse matrix
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

cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# STEP 2: Quality Control
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 3: Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# STEP 4: Feature Selection (Highly Variable Genes)
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

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# STEP 7: Clustering
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# STEP 8: UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP with cluster labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP - Seurat Clusters")

# ============================================================================
# STEP 9: Find Cluster Markers
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# ============================================================================
# STEP 10: Cell Type Annotation Based on Canonical Markers
# ============================================================================

# Define canonical marker genes for PBMC cell types
canonical_markers <- list(
  "CD4 T cells" = c("IL7R", "CCR7"),
  "CD14 Monocytes" = c("CD14", "LYZ"),
  "B cells" = c("MS4A1"),
  "CD8 T cells" = c("CD8A"),
  "NK cells" = c("GNLY", "NKG7")
)

# Visualize canonical markers on UMAP
for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  available_markers <- markers[markers %in% rownames(pbmc)]
  
  if (length(available_markers) > 0) {
    print(FeaturePlot(pbmc, features = available_markers, ncol = length(available_markers)) +
            plot_annotation(title = cell_type))
  }
}

# Visualize all canonical markers as dot plot
all_markers <- unique(unlist(canonical_markers))
available_all_markers <- all_markers[all_markers %in% rownames(pbmc)]

DotPlot(pbmc, features = available_all_markers) + 
  RotatedAxis() +
  ggtitle("Canonical PBMC Markers by Cluster")

# ============================================================================
# STEP 11: Automated Cell Type Assignment
# ============================================================================

# Calculate average expression of marker genes per cluster
cluster_avg <- AverageExpression(pbmc, features = available_all_markers, return.seurat = FALSE)
cluster_avg_matrix <- cluster_avg$RNA

# Function to assign cell types based on marker expression
assign_cell_type <- function(cluster_id, avg_expr) {
  scores <- numeric(length(canonical_markers))
  names(scores) <- names(canonical_markers)
  
  for (i in seq_along(canonical_markers)) {
    cell_type <- names(canonical_markers)[i]
    markers <- canonical_markers[[cell_type]]
    available <- markers[markers %in% rownames(avg_expr)]
    
    if (length(available) > 0) {
      scores[cell_type] <- mean(avg_expr[available, cluster_id])
    }
  }
  
  return(names(which.max(scores)))
}

# Assign cell types to each cluster
cluster_ids <- colnames(cluster_avg_matrix)
cell_type_assignments <- sapply(cluster_ids, function(x) assign_cell_type(x, cluster_avg_matrix))

# Create a mapping dataframe
cell_type_map <- data.frame(
  cluster = cluster_ids,
  cell_type = cell_type_assignments,
  stringsAsFactors = FALSE
)

print("Cell Type Assignments:")
print(cell_type_map)

# Add cell type annotations to Seurat object
current_clusters <- as.character(Idents(pbmc))
pbmc$cell_type <- cell_type_map$cell_type[match(current_clusters, cell_type_map$cluster)]

# ============================================================================
# STEP 12: Final Visualization with Cell Type Labels
# ============================================================================

# Set cell type as active identity
Idents(pbmc) <- "cell_type"

# UMAP with cell type labels
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5) + 
  ggtitle("UMAP -
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
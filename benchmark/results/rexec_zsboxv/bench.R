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
library(ggplot2)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

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

# ============================================================================
# 2. CREATE SEURAT OBJECT
# ============================================================================

# Create Seurat object with basic QC filtering
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# ============================================================================
# 3. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ============================================================================
# 4. NORMALIZATION
# ============================================================================

# Normalize the data using log-normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 5. FEATURE SELECTION
# ============================================================================

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Visualize top variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# 6. SCALING
# ============================================================================

# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 7. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 8. CLUSTERING
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ============================================================================
# 9. UMAP DIMENSIONALITY REDUCTION
# ============================================================================

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters on UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# 10. VISUALIZE CANONICAL PBMC MARKER GENES
# ============================================================================

# Define canonical PBMC marker genes
marker_genes <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                  "FCGR3A", "LYZ", "PPBP", "CD8A")

# Create FeaturePlot for all marker genes
FeaturePlot(pbmc, 
            features = marker_genes, 
            ncol = 3,
            pt.size = 0.5,
            order = TRUE)

# Optional: Create individual plots for each marker
for (gene in marker_genes) {
  p <- FeaturePlot(pbmc, features = gene, pt.size = 0.5, order = TRUE) +
    ggtitle(paste0(gene, " Expression"))
  print(p)
}

# Optional: Violin plots for marker genes by cluster
VlnPlot(pbmc, features = marker_genes, ncol = 3, pt.size = 0.1)

# Optional: Dot plot for marker genes
DotPlot(pbmc, features = marker_genes) + RotatedAxis()

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc2700_seurat.rds")

# Save UMAP plot
ggsave("umap_clusters.pdf", DimPlot(pbmc, reduction = "umap", label = TRUE), 
       width = 8, height = 6)

# Save marker gene feature plot
ggsave("marker_genes_featureplot.pdf", 
       FeaturePlot(pbmc, features = marker_genes, ncol = 3), 
       width = 12, height = 12)

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
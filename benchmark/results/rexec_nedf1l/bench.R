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
# 1. Load the data
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

# ============================================
# 2. Create Seurat object
# ============================================

# Create Seurat object with basic filtering
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3,      # Keep genes expressed in >= 3 cells
                           min.features = 200) # Keep cells with >= 200 genes

# ============================================
# 3. Quality control and filtering
# ============================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ============================================
# 4. Normalization and scaling
# ============================================

# Normalize the data (log-normalization)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data (z-score transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================
# 5. PCA (Principal Component Analysis)
# ============================================

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine number of PCs
ElbowPlot(pbmc, ndims = 50)

# ============================================
# 6. Build KNN graph and clustering
# ============================================

# Build KNN graph using 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Cluster cells with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check number of clusters
table(Idents(pbmc))

# ============================================
# 7. UMAP visualization
# ============================================

# Run UMAP using 10 PCs
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters on UMAP
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("PBMC Clusters (Resolution 0.5)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# Save the plot
ggsave("pbmc_umap_clusters.png", plot = p1, width = 10, height = 8, dpi = 300)

# ============================================
# 8. Additional visualizations
# ============================================

# UMAP colored by number of genes
p2 <- FeaturePlot(pbmc, features = "nFeature_RNA", reduction = "umap")
print(p2)

# UMAP colored by UMI counts
p3 <- FeaturePlot(pbmc, features = "nCount_RNA", reduction = "umap")
print(p3)

# ============================================
# 9. Summary statistics
# ============================================

cat("\n=== Analysis Summary ===\n")
cat("Total cells after filtering:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("\nCells per cluster:\n")
print(table(Idents(pbmc)))

# ============================================
# 10. Save the Seurat object
# ============================================

saveRDS(pbmc, file = "pbmc_analyzed.rds")
cat("\nSeurat object saved as 'pbmc_analyzed.rds'\n")

# To load later: pbmc <- readRDS("pbmc_analyzed.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
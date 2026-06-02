setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Complete Single-Cell RNA-seq Analysis Pipeline for PBMC Data
# This script performs end-to-end analysis of 2,700 PBMCs using Seurat

# Load required libraries
library(Seurat)
library(Matrix)
library(dplyr)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load Data
# ============================================================================
cat("Loading data...\n")

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

cat("Data loaded successfully!\n")
cat(sprintf("Matrix dimensions: %d genes x %d cells\n", nrow(matrix), ncol(matrix)))

# ============================================================================
# STEP 2: Create Seurat Object
# ============================================================================
cat("\nCreating Seurat object...\n")

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Keep genes detected in at least 3 cells
  min.features = 200  # Keep cells with at least 200 detected genes
)

cat(sprintf("Seurat object created with %d cells\n", ncol(pbmc)))

# ============================================================================
# STEP 3: Quality Control Metrics
# ============================================================================
cat("\nCalculating QC metrics...\n")

# Calculate percentage of mitochondrial genes
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics before filtering
cat("\nQC metrics before filtering:\n")
cat(sprintf("Median genes per cell: %.0f\n", median(pbmc$nFeature_RNA)))
cat(sprintf("Median UMIs per cell: %.0f\n", median(pbmc$nCount_RNA)))
cat(sprintf("Median percent MT: %.2f%%\n", median(pbmc$percent.mt)))

# ============================================================================
# STEP 4: Filter Low-Quality Cells
# ============================================================================
cat("\nFiltering low-quality cells...\n")
cat("Criteria: 200-2500 genes detected, <5% mitochondrial reads\n")

# Apply quality control filters
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & 
           nFeature_RNA < 2500 & 
           percent.mt < 5
)

cat(sprintf("Cells remaining after filtering: %d\n", ncol(pbmc)))

# ============================================================================
# STEP 5: Normalization
# ============================================================================
cat("\nNormalizing data...\n")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("Normalization complete!\n")

# ============================================================================
# STEP 6: Find Variable Features
# ============================================================================
cat("\nIdentifying highly variable features...\n")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

cat(sprintf("Identified %d variable features\n", length(VariableFeatures(pbmc))))

# ============================================================================
# STEP 7: Scale Data
# ============================================================================
cat("\nScaling data...\n")

# Get all genes for scaling
all.genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all.genes
)

cat("Data scaling complete!\n")

# ============================================================================
# STEP 8: Principal Component Analysis (PCA)
# ============================================================================
cat("\nRunning PCA...\n")

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = 50,
  verbose = FALSE
)

cat("PCA complete!\n")
cat(sprintf("Computed %d principal components\n", ncol(pbmc@reductions$pca)))

# ============================================================================
# STEP 9: Clustering
# ============================================================================
cat("\nPerforming clustering...\n")

# Find neighbors using first 10 PCs
pbmc <- FindNeighbors(
  pbmc,
  dims = 1:10
)

# Find clusters with resolution 0.5
pbmc <- FindClusters(
  pbmc,
  resolution = 0.5,
  verbose = FALSE
)

cat(sprintf("Identified %d clusters\n", length(unique(Idents(pbmc)))))
cat("Cluster sizes:\n")
print(table(Idents(pbmc)))

# ============================================================================
# STEP 10: UMAP Dimensionality Reduction
# ============================================================================
cat("\nComputing UMAP...\n")

pbmc <- RunUMAP(
  pbmc,
  dims = 1:10,
  verbose = FALSE
)

cat("UMAP complete!\n")

# ============================================================================
# STEP 11: Save Final Object
# ============================================================================
cat("\nSaving final Seurat object...\n")

saveRDS(pbmc, file = "pbmc_final.rds")

cat("Analysis complete! Seurat object saved as 'pbmc_final.rds'\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("Final number of cells: %d\n", ncol(pbmc)))
cat(sprintf("Final number of genes: %d\n", nrow(pbmc)))
cat(sprintf("Number of clusters: %d\n", length(unique(Idents(pbmc)))))
cat(sprintf("Variable features: %d\n", length(VariableFeatures(pbmc))))
cat(sprintf("PCA dimensions: %d\n", ncol(pbmc@reductions$pca)))
cat("\nCluster distribution:\n")
print(table(Idents(pbmc)))
cat(rep("=", 70), "\n", sep = "")

# Optional: Quick visualization (uncomment to generate plots)
# pdf("pbmc_analysis_plots.pdf", width = 10, height = 8)
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
# dev.off()

cat("\nTo load the object later, use: pbmc <- readRDS('pbmc_final.rds')\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
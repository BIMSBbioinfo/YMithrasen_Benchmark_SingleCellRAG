setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Complete Single-Cell RNA-seq Analysis Pipeline for PBMC Data
# This script performs end-to-end analysis of 10X Genomics scRNA-seq data

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data
# ============================================================================
cat("Loading data...\n")

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (feature identifiers)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(genes) <- c("gene_id", "gene_name")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_name
colnames(matrix) <- barcodes$barcode

cat("Data loaded successfully!\n")
cat(sprintf("Matrix dimensions: %d genes x %d cells\n", nrow(matrix), ncol(matrix)))

# ============================================================================
# STEP 2: Create Seurat object
# ============================================================================
cat("\nCreating Seurat object...\n")

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Keep genes detected in at least 3 cells
  min.features = 200  # Keep cells with at least 200 detected genes
)

cat(sprintf("Seurat object created with %d cells\n", ncol(pbmc)))

# ============================================================================
# STEP 3: Calculate QC metrics
# ============================================================================
cat("\nCalculating QC metrics...\n")

# Calculate percentage of mitochondrial genes
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Display QC metrics summary
cat("\nQC Metrics Summary:\n")
print(summary(pbmc$nFeature_RNA))
cat("\nMitochondrial percentage:\n")
print(summary(pbmc$percent.mt))

# ============================================================================
# STEP 4: Filter cells based on QC metrics
# ============================================================================
cat("\nFiltering cells...\n")
cat(sprintf("Cells before filtering: %d\n", ncol(pbmc)))

pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & 
           nFeature_RNA < 2500 & 
           percent.mt < 5
)

cat(sprintf("Cells after filtering: %d\n", ncol(pbmc)))

# ============================================================================
# STEP 5: Normalize the data
# ============================================================================
cat("\nNormalizing data...\n")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ============================================================================
# STEP 6: Find variable features
# ============================================================================
cat("\nFinding variable features...\n")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:\n")
print(top10)

# ============================================================================
# STEP 7: Scale the data
# ============================================================================
cat("\nScaling data...\n")

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 8: Perform PCA
# ============================================================================
cat("\nRunning PCA...\n")

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = 50,
  verbose = FALSE
)

cat("PCA completed. Top genes associated with PC1:\n")
print(pbmc[["pca"]], dims = 1, nfeatures = 5)

# ============================================================================
# STEP 9: Find neighbors and clusters
# ============================================================================
cat("\nFinding neighbors...\n")

pbmc <- FindNeighbors(
  pbmc,
  dims = 1:10,
  verbose = FALSE
)

cat("Clustering cells...\n")

pbmc <- FindClusters(
  pbmc,
  resolution = 0.5,
  verbose = FALSE
)

cat(sprintf("Number of clusters identified: %d\n", 
            length(unique(Idents(pbmc)))))

# ============================================================================
# STEP 10: Run UMAP
# ============================================================================
cat("\nRunning UMAP...\n")

pbmc <- RunUMAP(
  pbmc,
  dims = 1:10,
  verbose = FALSE
)

# ============================================================================
# STEP 11: Save the final Seurat object
# ============================================================================
cat("\nSaving Seurat object...\n")

saveRDS(pbmc, file = "pbmc_final.rds")

cat("\n============================================================================\n")
cat("Analysis complete!\n")
cat("============================================================================\n")
cat(sprintf("Final object contains %d cells and %d genes\n", 
            ncol(pbmc), nrow(pbmc)))
cat(sprintf("Number of clusters: %d\n", length(unique(Idents(pbmc)))))
cat("Saved to: pbmc_final.rds\n")
cat("\nTo load the object in a new session, use:\n")
cat("pbmc <- readRDS('pbmc_final.rds')\n")
cat("\nTo visualize clusters:\n")
cat("DimPlot(pbmc, reduction = 'umap', label = TRUE)\n")
cat("============================================================================\n")

# Optional: Print session info for reproducibility
cat("\nSession Info:\n")
sessionInfo()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
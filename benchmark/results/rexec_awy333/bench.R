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

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read gene information
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

cat("Initial dataset dimensions:\n")
cat("Cells:", ncol(pbmc), "\n")
cat("Genes:", nrow(pbmc), "\n\n")

# ============================================================================
# 2. QUALITY CONTROL METRICS
# ============================================================================

# Calculate mitochondrial percentage
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View QC metrics
cat("QC Metrics Summary (before filtering):\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Visualize QC metrics (optional)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)

# ============================================================================
# 3. QC FILTERING
# ============================================================================

# Apply filtering criteria:
# - nFeature_RNA > 200 and < 2500 (number of genes detected per cell)
# - percent.mt < 5 (mitochondrial content)

pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)

cat("\nDataset dimensions after QC filtering:\n")
cat("Cells:", ncol(pbmc), "\n")
cat("Genes:", nrow(pbmc), "\n\n")

# ============================================================================
# 4. NORMALIZATION
# ============================================================================

# Log-normalization
# Normalizes feature expression for each cell by total expression
# Multiplies by scale factor (10,000 by default)
# Log-transforms the result

pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

cat("Data normalized using log-normalization.\n\n")

# ============================================================================
# 5. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================

# Find the top 2000 most variable features
# These features exhibit high cell-to-cell variation

pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 most variable genes:\n")
print(top10)
cat("\n")

# Visualize variable features (optional)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# 6. SCALE DATA
# ============================================================================

# Scale the data
# - Shifts expression of each gene to mean = 0
# - Scales expression of each gene to variance = 1
# This prevents highly-expressed genes from dominating downstream analyses

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaled (mean = 0, variance = 1).\n\n")

# ============================================================================
# 7. SUMMARY AND SAVE
# ============================================================================

cat("Processing complete!\n")
cat("Final object contains:\n")
cat("- Cells:", ncol(pbmc), "\n")
cat("- Genes:", nrow(pbmc), "\n")
cat("- Variable features:", length(VariableFeatures(pbmc)), "\n\n")

# Save the processed Seurat object
saveRDS(pbmc, file = "pbmc_processed.rds")
cat("Processed data saved to: pbmc_processed.rds\n")

# Optional: View the Seurat object structure
print(pbmc)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
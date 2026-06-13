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
library(dplyr)

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (feature identifiers)
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
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics before filtering
cat("QC metrics before filtering:\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Optional: Create QC violin plots
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)

# Optional: Feature-feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# ============================================================================
# 3. QC FILTERING
# ============================================================================

# Apply QC filters
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

cat("\nDataset dimensions after QC filtering:\n")
cat("Cells:", ncol(pbmc), "\n")
cat("Genes:", nrow(pbmc), "\n")
cat("Cells removed:", ncol(matrix) - ncol(pbmc), "\n\n")

# ============================================================================
# 4. NORMALIZATION
# ============================================================================

# Log-normalization (default: scale factor = 10,000)
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

cat("Data normalized using log-normalization\n\n")

# ============================================================================
# 5. IDENTIFY VARIABLE FEATURES
# ============================================================================

# Find top 2000 highly variable features
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:\n")
print(top10)
cat("\n")

# Optional: Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 6. SCALE DATA
# ============================================================================

# Scale data (z-score transformation)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaled (z-score transformation applied)\n\n")

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("=== PROCESSING COMPLETE ===\n")
cat("Final Seurat object summary:\n")
print(pbmc)

# Save the processed Seurat object (optional)
# saveRDS(pbmc, file = "pbmc2700_processed.rds")

# To load later:
# pbmc <- readRDS("pbmc2700_processed.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
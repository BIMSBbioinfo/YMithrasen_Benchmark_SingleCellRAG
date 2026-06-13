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
colnames(genes) <- c("gene_id", "gene_name")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_name
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
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics before filtering
cat("QC metrics before filtering:\n")
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$percent.mt))

# Create QC violin plots
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1)

# Create feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Apply QC filters
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

cat("\nDataset dimensions after QC filtering:\n")
cat("Cells:", ncol(pbmc), "\n")
cat("Genes:", nrow(pbmc), "\n")
cat("Cells removed:", ncol(matrix) - ncol(pbmc), "\n\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

# Log-normalization (default: scale.factor = 10000)
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

cat("Data normalized using LogNormalize method\n\n")

# ============================================================================
# 4. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================

# Find top 2000 variable features
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 most variable genes:\n")
print(top10)
cat("\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 5. SCALE DATA
# ============================================================================

# Scale all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

cat("Data scaling complete\n")
cat("Scaled features:", length(all.genes), "\n\n")

# ============================================================================
# 6. SUMMARY AND SAVE
# ============================================================================

# Print final summary
cat("=== PROCESSING COMPLETE ===\n")
cat("Final dataset:\n")
cat("  Cells:", ncol(pbmc), "\n")
cat("  Genes:", nrow(pbmc), "\n")
cat("  Variable features:", length(VariableFeatures(pbmc)), "\n")
cat("  Mitochondrial threshold: < 5%\n")
cat("  Feature count range: 200-2500\n\n")

# Save the processed Seurat object
saveRDS(pbmc, file = "pbmc2700_processed.rds")
cat("Processed Seurat object saved as 'pbmc2700_processed.rds'\n")

# Optional: View the Seurat object structure
print(pbmc)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
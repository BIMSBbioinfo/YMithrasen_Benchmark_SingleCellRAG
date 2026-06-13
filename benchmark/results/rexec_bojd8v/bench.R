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

# Set working directory to where your files are located
# setwd("path/to/your/data")

# Load the Matrix Market format data
# Read the sparse matrix
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

# Calculate mitochondrial percentage
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View QC metrics
head(pbmc@meta.data)

# Summary statistics
cat("\nQC Metrics Summary:\n")
cat("Number of cells:", ncol(pbmc), "\n")
cat("Number of genes:", nrow(pbmc), "\n")
cat("\nnCount_RNA (UMI counts per cell):\n")
print(summary(pbmc$nCount_RNA))
cat("\nnFeature_RNA (genes detected per cell):\n")
print(summary(pbmc$nFeature_RNA))
cat("\npercent.mt (mitochondrial percentage):\n")
print(summary(pbmc$percent.mt))

# Create feature scatter plots for QC
# Plot 1: nCount_RNA vs percent.mt
plot1 <- FeatureScatter(pbmc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt") +
  ggtitle("QC: UMI Counts vs Mitochondrial Percentage") +
  theme_minimal()

# Plot 2: nCount_RNA vs nFeature_RNA
plot2 <- FeatureScatter(pbmc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA") +
  ggtitle("QC: UMI Counts vs Genes Detected") +
  theme_minimal()

# Display plots
print(plot1)
print(plot2)

# Optional: Save plots
ggsave("QC_nCount_vs_percentMT.png", plot1, width = 8, height = 6, dpi = 300)
ggsave("QC_nCount_vs_nFeature.png", plot2, width = 8, height = 6, dpi = 300)

# Optional: Create violin plots for additional QC visualization
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,
        pt.size = 0.1)

# Suggested filtering thresholds (adjust based on your data)
cat("\n\nSuggested filtering criteria:\n")
cat("- Filter cells with nFeature_RNA < 200 or > 2500\n")
cat("- Filter cells with percent.mt > 5%\n")
cat("\nTo apply filtering, use:\n")
cat("pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
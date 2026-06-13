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

# Load the data files
# Read the sparse matrix in Matrix Market format
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

# Calculate mitochondrial percentage
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View basic statistics
print(paste("Number of cells:", ncol(pbmc)))
print(paste("Number of genes:", nrow(pbmc)))
print(head(pbmc@meta.data))

# Summary statistics
summary(pbmc$nCount_RNA)
summary(pbmc$nFeature_RNA)
summary(pbmc$percent.mt)

# Create quality control plots
# Plot 1: nCount_RNA vs percent.mt
plot1 <- FeatureScatter(pbmc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt") +
  ggtitle("UMI Counts vs Mitochondrial Percentage") +
  theme_minimal()

# Plot 2: nCount_RNA vs nFeature_RNA
plot2 <- FeatureScatter(pbmc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA") +
  ggtitle("UMI Counts vs Number of Genes") +
  theme_minimal()

# Display plots
print(plot1)
print(plot2)

# Optional: Save plots
ggsave("QC_nCount_vs_percent_mt.png", plot1, width = 8, height = 6, dpi = 300)
ggsave("QC_nCount_vs_nFeature.png", plot2, width = 8, height = 6, dpi = 300)

# Optional: Create violin plots for additional QC visualization
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Print QC metrics summary
cat("\n=== Quality Control Metrics ===\n")
cat("Total UMI counts per cell:\n")
print(summary(pbmc$nCount_RNA))
cat("\nNumber of genes detected per cell:\n")
print(summary(pbmc$nFeature_RNA))
cat("\nMitochondrial percentage per cell:\n")
print(summary(pbmc$percent.mt))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
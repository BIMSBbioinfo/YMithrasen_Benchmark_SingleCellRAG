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

# Load the sparse matrix data
matrix <- readMM("matrix.mtx")
genes <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes[, 2]  # Use HGNC gene symbols
colnames(matrix) <- barcodes[, 1]  # Use cell barcodes

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
print(summary(pbmc$nCount_RNA))
print(summary(pbmc$nFeature_RNA))
print(summary(pbmc$percent.mt))

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
ggsave("QC_nCount_vs_percentMT.png", plot1, width = 8, height = 6, dpi = 300)
ggsave("QC_nCount_vs_nFeature.png", plot2, width = 8, height = 6, dpi = 300)

# Optional: Create violin plots for additional QC visualization
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Print QC metrics summary
cat("\n=== Quality Control Metrics ===\n")
cat("Median UMI counts per cell:", median(pbmc$nCount_RNA), "\n")
cat("Median genes per cell:", median(pbmc$nFeature_RNA), "\n")
cat("Median mitochondrial percentage:", median(pbmc$percent.mt), "%\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
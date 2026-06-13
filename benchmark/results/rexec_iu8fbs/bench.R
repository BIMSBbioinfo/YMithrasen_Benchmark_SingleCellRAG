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

# Read the 10X format data files
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

# Calculate percentage of mitochondrial genes
# Mitochondrial genes start with "MT-" in human data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View the metadata to confirm calculations
head(pbmc@meta.data)

# Create violin plots for QC metrics
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,
        pt.size = 0.1)

# Optional: Save the plot
ggsave("QC_violin_plot.png", width = 12, height = 4, dpi = 300)

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\nNumber of cells:", ncol(pbmc), "\n")
cat("Number of genes:", nrow(pbmc), "\n\n")

cat("nFeature_RNA (genes per cell):\n")
print(summary(pbmc$nFeature_RNA))

cat("\nnCount_RNA (UMI counts per cell):\n")
print(summary(pbmc$nCount_RNA))

cat("\npercent.mt (mitochondrial percentage):\n")
print(summary(pbmc$percent.mt))

# Optional: Create additional scatter plots to visualize relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Display scatter plots side by side
plot1 + plot2

# Optional: Save scatter plots
ggsave("QC_scatter_plots.png", plot = plot1 + plot2, width = 12, height = 5, dpi = 300)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
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
library(dplyr)

# Set working directory (adjust to your file location)
# setwd("path/to/your/data")

# ============================================
# 1. Load the 10x data files
# ============================================

# Read the sparse matrix (UMI counts)
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

# ============================================
# 2. Create Seurat object
# ============================================

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# ============================================
# 3. Calculate mitochondrial gene percentage
# ============================================

# Identify mitochondrial genes (genes starting with "MT-")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Alternative manual calculation if needed:
# mito_genes <- grep("^MT-", rownames(pbmc), value = TRUE)
# pbmc[["percent.mt"]] <- colSums(pbmc@assays$RNA@counts[mito_genes, ]) / 
#                          colSums(pbmc@assays$RNA@counts) * 100

# ============================================
# 4. View QC metrics summary
# ============================================

# Display summary statistics
cat("QC Metrics Summary:\n")
cat("==================\n")
summary(pbmc$nFeature_RNA)
summary(pbmc$nCount_RNA)
summary(pbmc$percent.mt)

# View first few cells
head(pbmc@meta.data)

# ============================================
# 5. Visualize QC metrics with violin plots
# ============================================

# Create violin plot for all three metrics
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1)

# Save the plot
ggsave("QC_violin_plots.png", width = 12, height = 4, dpi = 300)

# ============================================
# 6. Additional visualizations (optional)
# ============================================

# Scatter plots to visualize relationships between metrics
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Combine plots
library(patchwork)
plot1 + plot2

# Save combined scatter plots
ggsave("QC_scatter_plots.png", width = 10, height = 4, dpi = 300)

# ============================================
# 7. Identify potential filtering thresholds
# ============================================

cat("\nSuggested QC filtering thresholds:\n")
cat("===================================\n")
cat("Low quality cells (low nFeature_RNA): < 200\n")
cat("Potential doublets (high nFeature_RNA): >", quantile(pbmc$nFeature_RNA, 0.98), "\n")
cat("High mitochondrial content: >", quantile(pbmc$percent.mt, 0.95), "%\n")

# ============================================
# 8. Apply filtering (example thresholds)
# ============================================

# Filter cells based on QC metrics
pbmc_filtered <- subset(pbmc, 
                        subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 2500 & 
                                percent.mt < 5)

cat("\nCells before filtering:", ncol(pbmc), "\n")
cat("Cells after filtering:", ncol(pbmc_filtered), "\n")
cat("Cells removed:", ncol(pbmc) - ncol(pbmc_filtered), "\n")

# Visualize filtered data
VlnPlot(pbmc_filtered, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0.1)

ggsave("QC_violin_plots_filtered.png", width = 12, height = 4, dpi = 300)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
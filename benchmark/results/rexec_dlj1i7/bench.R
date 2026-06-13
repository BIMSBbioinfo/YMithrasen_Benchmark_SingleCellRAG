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

# Read the data files
# Read the sparse matrix in Matrix Market format
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

# Create a Seurat object
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
cat("\nnFeature_RNA (genes per cell):\n")
summary(pbmc$nFeature_RNA)

cat("\nnCount_RNA (UMI counts per cell):\n")
summary(pbmc$nCount_RNA)

cat("\npercent.mt (mitochondrial percentage):\n")
summary(pbmc$percent.mt)

# Optional: Create a more detailed visualization with ggplot2
library(reshape2)

# Extract metadata
qc_data <- pbmc@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]

# Melt data for ggplot2
qc_melt <- melt(qc_data)

# Create custom violin plot
ggplot(qc_melt, aes(x = variable, y = value, fill = variable)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_jitter(alpha = 0.1, size = 0.5) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Quality Control Metrics for 2,700 PBMCs",
       x = "Metric",
       y = "Value")

ggsave("QC_violin_plot_detailed.png", width = 12, height = 4, dpi = 300)

# Optional: Filter cells based on QC metrics (common thresholds)
# Uncomment and adjust thresholds as needed
# pbmc_filtered <- subset(pbmc, 
#                         subset = nFeature_RNA > 200 & 
#                                  nFeature_RNA < 2500 & 
#                                  percent.mt < 5)
# 
# cat("\nCells before filtering:", ncol(pbmc))
# cat("\nCells after filtering:", ncol(pbmc_filtered))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
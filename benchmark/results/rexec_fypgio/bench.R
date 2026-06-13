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
# 1. Load the data files
# ============================================

# Read the sparse matrix (Matrix Market format)
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================
# 2. Assign row and column names to matrix
# ============================================

# Set gene symbols as row names and barcodes as column names
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# ============================================
# 3. Create Seurat object
# ============================================

pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# ============================================
# 4. Calculate mitochondrial gene percentage
# ============================================

# Identify mitochondrial genes (genes starting with "MT-")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Alternative manual calculation if needed:
# mito_genes <- grep("^MT-", rownames(pbmc), value = TRUE)
# pbmc[["percent.mt"]] <- colSums(pbmc@assays$RNA@counts[mito_genes, ]) / 
#                          colSums(pbmc@assays$RNA@counts) * 100

# ============================================
# 5. View QC metrics summary
# ============================================

# Display summary statistics
cat("\n=== QC Metrics Summary ===\n")
summary(pbmc$nFeature_RNA)
summary(pbmc$nCount_RNA)
summary(pbmc$percent.mt)

# View first few cells
head(pbmc@meta.data)

# ============================================
# 6. Create violin plots for QC metrics
# ============================================

# Violin plot using Seurat's built-in function
p1 <- VlnPlot(
  pbmc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1  # Point size for individual cells
)

# Display the plot
print(p1)

# Save the plot
ggsave("QC_violin_plots.png", plot = p1, width = 12, height = 4, dpi = 300)

# ============================================
# 7. Alternative: Custom ggplot2 violin plots
# ============================================

# Extract metadata for custom plotting
metadata <- pbmc@meta.data

# Create individual plots with more customization
p2 <- ggplot(metadata, aes(x = "", y = nFeature_RNA)) +
  geom_violin(fill = "#C77CFF", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.3) +
  labs(title = "Number of Genes per Cell", 
       x = "", 
       y = "nFeature_RNA") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(metadata, aes(x = "", y = nCount_RNA)) +
  geom_violin(fill = "#00BFC4", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.3) +
  labs(title = "Total UMI Counts per Cell", 
       x = "", 
       y = "nCount_RNA") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(metadata, aes(x = "", y = percent.mt)) +
  geom_violin(fill = "#F8766D", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.3) +
  labs(title = "Mitochondrial Gene Percentage", 
       x = "", 
       y = "% Mitochondrial Genes") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine plots using patchwork (optional)
# library(patchwork)
# combined_plot <- p2 | p3 | p4
# ggsave("QC_custom_violin_plots.png", plot = combined_plot, width = 12, height = 4, dpi = 300)

# ============================================
# 8. Additional QC visualizations
# ============================================

# Scatter plots to identify outliers
p5 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p6 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

print(p5)
print(p6)

# ============================================
# 9. Suggested filtering thresholds (adjust based on your data)
# ============================================

cat("\n=== Suggested Filtering Criteria ===\n")
cat("Based on typical PBMC QC standards:\n")
cat("- nFeature_RNA: 200 - 2500 genes\n")
cat("- percent.mt: < 5%\n")
cat("\nCells before filtering:", ncol(pbmc), "\n")

# Apply filters (uncomment to execute)
# pbmc_filtered <- subset(pbmc, 
#                         subset = nFeature_RNA > 200 & 
#                                  nFeature_RNA < 2500 & 
#                                  percent.mt < 5)
# 
# cat("Cells after filtering:", ncol(pbmc_filtered), "\n")

# ============================================
# 10. Save the Seurat object
# ============================================

saveRDS(pbmc, file = "pbmc_2700_qc.rds")
cat("\nSeurat object saved as 'pbmc_2700_qc.rds'\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
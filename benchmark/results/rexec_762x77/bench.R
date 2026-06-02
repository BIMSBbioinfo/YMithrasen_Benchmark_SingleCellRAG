setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
Here's executable R code to perform quality control on your 10x Genomics single-cell RNA-seq data:

```r
# Load required libraries
library(Matrix)
library(ggplot2)
library(dplyr)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ===== 1. Load the data =====
cat("Loading data...\n")

# Read the sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell IDs)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

cat("Data loaded successfully!\n")
cat("Dimensions:", nrow(matrix), "genes x", ncol(matrix), "cells\n\n")

# ===== 2. Calculate QC metrics =====
cat("Calculating QC metrics...\n")

# Number of genes detected per cell (genes with UMI > 0)
n_genes_per_cell <- colSums(matrix > 0)

# Total UMI counts per cell
n_umis_per_cell <- colSums(matrix)

# Identify mitochondrial genes (genes starting with "MT-" or "mt-")
mito_genes <- grep("^MT-|^mt-", rownames(matrix), value = TRUE)
cat("Found", length(mito_genes), "mitochondrial genes\n")

# Calculate percentage of mitochondrial UMIs per cell
if (length(mito_genes) > 0) {
  mito_counts <- colSums(matrix[mito_genes, , drop = FALSE])
  percent_mito <- (mito_counts / n_umis_per_cell) * 100
} else {
  percent_mito <- rep(0, ncol(matrix))
  cat("Warning: No mitochondrial genes found!\n")
}

# Create QC dataframe
qc_metrics <- data.frame(
  barcode = barcodes$barcode,
  n_genes = n_genes_per_cell,
  n_umis = n_umis_per_cell,
  percent_mito = percent_mito
)

# ===== 3. Visualize QC metrics =====
cat("\nGenerating QC plots...\n")

# Plot 1: Number of genes per cell
p1 <- ggplot(qc_metrics, aes(x = n_genes)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 200, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 2500, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Genes Detected per Cell", x = "Number of Genes", y = "Number of Cells") +
  theme_minimal()

# Plot 2: Total UMI counts per cell
p2 <- ggplot(qc_metrics, aes(x = n_umis)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "black") +
  geom_vline(xintercept = 500, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Total UMI Counts per Cell", x = "UMI Counts", y = "Number of Cells") +
  theme_minimal()

# Plot 3: Percentage of mitochondrial genes
p3 <- ggplot(qc_metrics, aes(x = percent_mito)) +
  geom_histogram(bins = 50, fill = "coral", color = "black") +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Mitochondrial Gene Percentage", x = "% Mitochondrial", y = "Number of Cells") +
  theme_minimal()

# Plot 4: Scatter plot - genes vs UMIs
p4 <- ggplot(qc_metrics, aes(x = n_umis, y = n_genes, color = percent_mito)) +
  geom_point(alpha = 0.5) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Genes vs UMIs", x = "Total UMI Counts", y = "Number of Genes", 
       color = "% Mito") +
  theme_minimal()

# Save plots
ggsave("qc_genes_per_cell.png", p1, width = 8, height = 6)
ggsave("qc_umi_counts.png", p2, width = 8, height = 6)
ggsave("qc_mito_percentage.png", p3, width = 8, height = 6)
ggsave("qc_genes_vs_umis.png", p4, width = 8, height = 6)

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)

# ===== 4. Print summary statistics =====
cat("\n===== QC Summary Statistics =====\n")
cat("Total cells before filtering:", ncol(matrix), "\n\n")

cat("Number of genes per cell:\n")
print(summary(n_genes_per_cell))

cat("\nTotal UMI counts per cell:\n")
print(summary(n_umis_per_cell))

cat("\nPercentage of mitochondrial genes:\n")
print(summary(percent_mito))

# ===== 5. Filter low-quality cells =====
cat("\n===== Filtering Cells =====\n")

# Define filtering thresholds
min_genes <- 200      # Minimum genes detected
max_genes <- 2500     # Maximum genes (to remove potential doublets)
min_umis <- 500       # Minimum UMI counts
max_mito <- 5         # Maximum percentage of mitochondrial genes

cat("Filtering criteria:\n")
cat("  - Min genes per cell:", min_genes, "\n")
cat("  - Max genes per cell:", max_genes, "\n")
cat("  - Min UMI counts:", min_umis, "\n")
cat("  - Max mitochondrial %:", max_mito, "\n\n")

# Apply filters
cells_to_keep <- (qc_metrics$n_genes >= min_genes) & 
                 (qc_metrics$n_genes <= max_genes) & 
                 (qc_metrics$n_umis >= min_umis) & 
                 (qc_metrics$percent_mito <= max_mito)

# Filter the matrix
filtered_matrix <- matrix[, cells_to_keep]
filtered_qc <- qc_metrics[cells_to_keep, ]

cat("Cells passing QC:", sum(cells_to_keep), "/", ncol(matrix), 
    "(", round(100 * sum(cells_to_keep) / ncol(matrix), 2), "%)\n")
cat("Cells removed:", sum(!cells_to_keep), "\n\n")

# ===== 6. Additional gene filtering =====
cat("===== Filtering Genes =====\n")

# Remove genes not detected in any cell
genes_detected <- rowSums(filtered_matrix > 0)
genes_to_keep <- genes_detected >= 3  # Keep genes detected in at least 3 cells

filtered_matrix <- filtered_matrix[genes_to_keep, ]

cat("Genes passing filter:", sum(genes_to_keep), "/", nrow(matrix), "\n")
cat("Final matrix dimensions:", nrow(filtered_matrix), "genes x", 
    ncol(filtered_matrix), "cells\n\n")

# ===== 7. Save
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
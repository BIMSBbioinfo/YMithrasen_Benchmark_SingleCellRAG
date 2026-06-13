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

# Set working directory (adjust to your data location)
# setwd("path/to/your/data")

# ============================================================================
# Step 1: Load the data files
# ============================================================================

# Read the sparse matrix (Matrix Market format)
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================================================
# Step 2: Format the matrix with proper row and column names
# ============================================================================

# Assign gene symbols as row names and barcodes as column names
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix)[1], "genes x", dim(matrix)[2], "cells\n")

# ============================================================================
# Step 3: Create Seurat object
# ============================================================================

# Create Seurat object with basic QC filtering
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,      # Include genes detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected genes
)

# Display basic information
print(pbmc)

# ============================================================================
# Step 4: Basic QC metrics (optional but recommended)
# ============================================================================

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# View QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ============================================================================
# Step 5: Normalize the data
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# Step 6: Identify highly variable features using VST method
# ============================================================================

pbmc <- FindVariableFeatures(
  pbmc, 
  selection.method = "vst",  # Variance Stabilizing Transformation
  nfeatures = 2000           # Select top 2000 features
)

# Get the top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("\nTop 10 highly variable genes:\n")
print(top10)

# ============================================================================
# Step 7: Visualize variable features with top 10 highlighted
# ============================================================================

# Create the variable feature plot
plot1 <- VariableFeaturePlot(pbmc)

# Add labels for top 10 genes
plot2 <- LabelPoints(
  plot = plot1, 
  points = top10, 
  repel = TRUE,           # Use ggrepel to avoid label overlap
  xnudge = 0,
  ynudge = 0
)

# Display the plot
print(plot2)

# Save the plot
ggsave("variable_features_plot.png", plot = plot2, width = 10, height = 6, dpi = 300)

# ============================================================================
# Additional visualization: Alternative plot style
# ============================================================================

# You can also create a custom plot with more control
hvf_data <- HVFInfo(pbmc)
hvf_data$gene <- rownames(hvf_data)
hvf_data$variable <- rownames(hvf_data) %in% VariableFeatures(pbmc)
hvf_data$top10 <- rownames(hvf_data) %in% top10

# Custom ggplot
custom_plot <- ggplot(hvf_data, aes(x = mean, y = variance.standardized)) +
  geom_point(aes(color = variable), alpha = 0.5, size = 1) +
  geom_point(data = subset(hvf_data, top10), color = "red", size = 2) +
  geom_text(data = subset(hvf_data, top10), 
            aes(label = gene), 
            hjust = -0.1, vjust = 0.5, size = 3.5, color = "red") +
  scale_color_manual(values = c("grey", "blue")) +
  labs(
    title = "Highly Variable Features (VST method)",
    subtitle = paste("Top 2000 features selected, Top 10 labeled in red"),
    x = "Average Expression",
    y = "Standardized Variance",
    color = "Variable Feature"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(custom_plot)
ggsave("variable_features_custom.png", plot = custom_plot, width = 10, height = 6, dpi = 300)

# ============================================================================
# Summary statistics
# ============================================================================

cat("\n=== Summary ===\n")
cat("Total cells after filtering:", ncol(pbmc), "\n")
cat("Total genes after filtering:", nrow(pbmc), "\n")
cat("Number of variable features:", length(VariableFeatures(pbmc)), "\n")
cat("\nTop 10 variable genes:\n")
print(top10)

# Save the Seurat object for future use
saveRDS(pbmc, file = "pbmc_seurat_object.rds")
cat("\nSeurat object saved as 'pbmc_seurat_object.rds'\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
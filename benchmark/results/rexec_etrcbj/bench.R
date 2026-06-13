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

# Set working directory to where your files are located
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data files
# ============================================================================

# Read the sparse matrix in Matrix Market format
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# ============================================================================
# STEP 2: Format the matrix with proper row and column names
# ============================================================================

# Assign gene symbols as row names and barcodes as column names
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Check dimensions
cat("Matrix dimensions:", dim(matrix), "\n")
cat("Number of genes:", nrow(matrix), "\n")
cat("Number of cells:", ncol(matrix), "\n")

# ============================================================================
# STEP 3: Create Seurat object
# ============================================================================

# Create Seurat object with basic QC filtering
# min.cells: include genes detected in at least 3 cells
# min.features: include cells with at least 200 detected genes
pbmc <- CreateSeuratObject(
  counts = matrix,
  project = "PBMC_2700",
  min.cells = 3,
  min.features = 200
)

# View basic information
print(pbmc)

# ============================================================================
# STEP 4: Quality Control (Optional but recommended)
# ============================================================================

# Calculate percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics (adjust thresholds as needed)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("Cells after QC filtering:", ncol(pbmc), "\n")

# ============================================================================
# STEP 5: Normalize the data
# ============================================================================

# Log-normalization: normalize gene expression for each cell by total expression,
# multiply by scale factor (10,000), and log-transform
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("Data normalized successfully\n")

# ============================================================================
# STEP 6: Identify highly variable genes
# ============================================================================

# Find the top 2000 most variable genes
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", top10, "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# STEP 7: Scale the data
# ============================================================================

# Scale data so that mean expression = 0 and variance = 1 across cells
# This is done for all genes by default, but you can specify features
all_genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all_genes,
  verbose = TRUE
)

cat("Data scaled successfully\n")

# Alternative: Scale only variable features (faster, recommended for large datasets)
# pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# ============================================================================
# STEP 8: Verify the results
# ============================================================================

# Check the scaled data
scaled_data <- GetAssayData(pbmc, slot = "scale.data")
cat("Scaled data dimensions:", dim(scaled_data), "\n")

# Verify mean ≈ 0 and variance ≈ 1 for a sample gene
sample_gene <- VariableFeatures(pbmc)[1]
gene_values <- scaled_data[sample_gene, ]
cat("\nVerification for gene:", sample_gene, "\n")
cat("Mean:", mean(gene_values), "\n")
cat("Variance:", var(gene_values), "\n")

# ============================================================================
# STEP 9: Save the processed Seurat object (optional)
# ============================================================================

saveRDS(pbmc, file = "pbmc_processed.rds")
cat("\nProcessed Seurat object saved as 'pbmc_processed.rds'\n")

# To load later: pbmc <- readRDS("pbmc_processed.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
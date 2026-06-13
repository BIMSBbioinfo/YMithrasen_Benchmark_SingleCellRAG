setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(DoubletFinder)
library(Matrix)
library(ggplot2)
library(dplyr)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes and genes
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes$V2  # Use HGNC symbols
colnames(matrix) <- barcodes$V1

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Loaded data with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Use first 10 PCs (adjust based on ElbowPlot)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ============================================================================
# 8. UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)

# ============================================================================
# 9. DOUBLET DETECTION WITH DOUBLETFINDER
# ============================================================================

cat("\n=== Running DoubletFinder ===\n")

# Pre-process for DoubletFinder
# pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select optimal pK
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
cat("Optimal pK:", pK_optimal, "\n")

# Homotypic Doublet Proportion Estimate
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
cat("Homotypic doublet proportion:", homotypic.prop, "\n")

# Assuming 7.5% doublet formation rate (typical for ~8K cells loaded)
# Adjust based on your experimental setup
nExp_poi <- round(0.075 * ncol(pbmc))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

cat("Expected doublets:", nExp_poi, "\n")
cat("Adjusted expected doublets:", nExp_poi.adj, "\n")

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, 
                         PCs = 1:10, 
                         pN = 0.25, 
                         pK = pK_optimal, 
                         nExp = nExp_poi.adj, 
                         reuse.pANN = FALSE, 
                         sct = FALSE)

# Get the DoubletFinder classification column name (it's auto-generated)
df_column <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)
cat("DoubletFinder column:", df_column, "\n")

# Rename for easier access
pbmc@meta.data$doublet_classification <- pbmc@meta.data[[df_column]]

# ============================================================================
# 10. VISUALIZATION
# ============================================================================

# Count singlets and doublets
table(pbmc@meta.data$doublet_classification)

# UMAP colored by doublet classification
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "doublet_classification", 
              cols = c("Singlet" = "#3498db", "Doublet" = "#e74c3c")) +
  ggtitle("PBMC 2700: Singlets vs Doublets") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# UMAP with clusters and doublets side by side
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clusters") +
  theme_minimal()

p3 <- DimPlot(pbmc, reduction = "umap", group.by = "doublet_classification",
              cols = c("Singlet" = "#3498db", "Doublet" = "#e74c3c")) +
  ggtitle("Doublet Classification") +
  theme_minimal()

# Display plots
print(p1)
print(p2 | p3)  # Side by side using patchwork

# Save plots
ggsave("umap_doublets.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("umap_clusters_and_doublets.png", plot = p2 | p3, width = 14, height = 6, dpi = 300)

# ============================================================================
# 11. SUMMARY STATISTICS
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Total cells:", ncol(pbmc), "\n")
cat("Singlets:", sum(pbmc@meta.data$doublet_classification == "Singlet"), "\n")
cat("Doublets:", sum(pbmc@meta.data$doublet_classification == "Doublet"), "\n")
cat("Doublet rate:", 
    round(100 * sum(pbmc@meta.data$doublet_classification == "Doublet") / ncol(pbmc), 2), "%\n")

# Doublets per cluster
doublet_by_cluster <- table(pbmc@meta.data$seurat_clusters, 
                            pbmc@meta.data$doublet_classification)
print(doublet_by_cluster)

# ============================================================================
# 12. OPTIONAL: REMOVE DOUBLETS AND RE-ANALYZE
# ============================================================================

# Subset to keep only singlets
pbmc_singlets <- subset(pbmc, subset = doublet_classification == "Singlet")

cat("\nSinglet-only dataset:", ncol(pbmc_singlets), "cells\n")

# You can now re-run clustering on singlets if desired
# pbmc_singlets <- FindNeighbors(pbmc_singlets, dims = 1:10)
# pbmc_singlets <- FindClusters(pbmc_singlets, resolution = 0.5)
# pbmc_singlets <- RunUMAP(pbmc_singlets, dims = 1:10)

# Save the Seurat object
saveRDS(pbmc, file = "pbmc2700_with_doublets.rds")
saveRDS(pbmc_singlets, file = "pbmc2700_singlets_only.rds")

cat("\nAnalysis complete!\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
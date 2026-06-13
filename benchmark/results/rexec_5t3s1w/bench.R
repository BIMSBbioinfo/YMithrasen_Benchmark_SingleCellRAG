setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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

# Estimate doublet rate (assuming ~7.5% for 2700 cells)
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

# Find the DoubletFinder classification column
df_col <- grep("DF.classifications", colnames(pbmc@meta.data), value = TRUE)
pANN_col <- grep("pANN", colnames(pbmc@meta.data), value = TRUE)

# Rename for easier access
pbmc@meta.data$doublet_classification <- pbmc@meta.data[[df_col]]
pbmc@meta.data$doublet_score <- pbmc@meta.data[[pANN_col]]

# ============================================================================
# 10. VISUALIZATION
# ============================================================================

# Count singlets and doublets
doublet_summary <- table(pbmc@meta.data$doublet_classification)
cat("\nDoublet Detection Summary:\n")
print(doublet_summary)
cat("Doublet rate:", round(doublet_summary["Doublet"] / sum(doublet_summary) * 100, 2), "%\n")

# UMAP colored by doublet classification
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "doublet_classification", 
              cols = c("Singlet" = "#3498db", "Doublet" = "#e74c3c")) +
  ggtitle("PBMC 2700: Singlets vs Doublets") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# UMAP colored by doublet score
p2 <- FeaturePlot(pbmc, reduction = "umap", features = "doublet_score") +
  ggtitle("DoubletFinder Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# UMAP colored by cluster with doublets highlighted
p3 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Combined visualization
library(patchwork)
combined_plot <- (p1 | p2) / p3
print(combined_plot)

# Save the plot
ggsave("PBMC_doublet_analysis.png", combined_plot, width = 14, height = 10, dpi = 300)

# ============================================================================
# 11. OPTIONAL: REMOVE DOUBLETS AND RE-ANALYZE
# ============================================================================

# Create a singlet-only object
pbmc_singlets <- subset(pbmc, subset = doublet_classification == "Singlet")

cat("\nSinglets only:", ncol(pbmc_singlets), "cells\n")

# You can now re-run clustering on singlets if desired
# pbmc_singlets <- FindNeighbors(pbmc_singlets, dims = 1:10)
# pbmc_singlets <- FindClusters(pbmc_singlets, resolution = 0.5)
# pbmc_singlets <- RunUMAP(pbmc_singlets, dims = 1:10)
# DimPlot(pbmc_singlets, reduction = "umap", label = TRUE)

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_with_doublets.rds")
saveRDS(pbmc_singlets, file = "pbmc_singlets_only.rds")

# Export metadata
write.csv(pbmc@meta.data, "pbmc_metadata_with_doublets.csv", row.names = TRUE)

cat("\n=== Analysis Complete ===\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
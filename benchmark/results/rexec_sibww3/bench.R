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

print(paste("Loaded", ncol(pbmc), "cells and", nrow(pbmc), "genes"))

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

print(paste("After QC filtering:", ncol(pbmc), "cells remaining"))

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
print("Top 10 variable genes:")
print(top10)

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Use first 10 PCs (adjust based on elbow plot)
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

# Pre-process for DoubletFinder
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Visualize pK selection
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  ggtitle("pK Selection")

# Select optimal pK (highest BCmetric)
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
print(paste("Optimal pK:", pK_optimal))

# Homotypic Doublet Proportion Estimate
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# Expected doublet rate (adjust based on your experiment)
# For 2,700 cells, expected doublet rate is ~2-3%
nExp_poi <- round(0.023 * ncol(pbmc))  # 2.3% expected doublet rate
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

print(paste("Expected doublets:", nExp_poi))
print(paste("Adjusted for homotypic doublets:", nExp_poi.adj))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, 
                         PCs = 1:10, 
                         pN = 0.25, 
                         pK = pK_optimal, 
                         nExp = nExp_poi.adj, 
                         reuse.pANN = FALSE, 
                         sct = FALSE)

# Get the DoubletFinder column name (it's auto-generated)
df_column <- colnames(pbmc@meta.data)[grep("DF.classifications", colnames(pbmc@meta.data))]

# Rename for easier access
pbmc@meta.data$doublet_classification <- pbmc@meta.data[[df_column]]

# ============================================================================
# 10. VISUALIZATION
# ============================================================================

# UMAP colored by doublet classification
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "doublet_classification", 
              cols = c("Singlet" = "#3498db", "Doublet" = "#e74c3c")) +
  ggtitle("Doublet Detection Results") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# UMAP with clusters and doublets highlighted
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat Clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Side-by-side comparison
library(patchwork)
p1 + p2

# Summary statistics
doublet_summary <- table(pbmc@meta.data$doublet_classification)
print("Doublet Detection Summary:")
print(doublet_summary)
print(paste("Doublet rate:", 
            round(doublet_summary["Doublet"] / sum(doublet_summary) * 100, 2), "%"))

# Doublets per cluster
doublet_by_cluster <- table(pbmc@meta.data$seurat_clusters, 
                            pbmc@meta.data$doublet_classification)
print("Doublets by cluster:")
print(doublet_by_cluster)

# ============================================================================
# 11. OPTIONAL: REMOVE DOUBLETS AND RE-ANALYZE
# ============================================================================

# Create a singlet-only object
pbmc_singlets <- subset(pbmc, subset = doublet_classification == "Singlet")

print(paste("Singlets retained:", ncol(pbmc_singlets), "cells"))

# Visualize singlets only
DimPlot(pbmc_singlets, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP - Singlets Only") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_with_doublets.rds")
saveRDS(pbmc_singlets, file = "pbmc_singlets_only.rds")

# Export metadata
write.csv(pbmc@meta.data, file = "pbmc_metadata_with_doublets.csv")

# Save plots
ggsave("umap_doublets.pdf", plot = p1, width = 8, height = 6)
ggsave("umap_clusters.pdf", plot = p2, width = 8, height = 6)

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
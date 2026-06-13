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
library(dplyr)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format sparse matrix
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

cat("Loaded data with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# 2. QUALITY CONTROL AND PREPROCESSING
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After filtering:", ncol(pbmc), "cells remaining\n")

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
# 6. PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)

# Visualize PCA
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 8. UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)

# ============================================================================
# 9. CELL CYCLE SCORING
# ============================================================================

# Load cell cycle genes (built into Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Perform cell cycle scoring
pbmc <- CellCycleScoring(pbmc, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = TRUE)

# Check the cell cycle phases assigned
table(pbmc$Phase)

# ============================================================================
# 10. VISUALIZE CELL CYCLE PHASES ON UMAP
# ============================================================================

# UMAP colored by cell cycle phase
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Phase", 
              cols = c("G1" = "#F8766D", "G2M" = "#00BA38", "S" = "#619CFF")) +
  ggtitle("PBMC Cells Colored by Cell Cycle Phase") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# UMAP colored by clusters with cell cycle phase overlay
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("PBMC Clusters")

print(p2)

# Feature plots for cell cycle scores
p3 <- FeaturePlot(pbmc, features = c("S.Score", "G2M.Score"), 
                  blend = FALSE, ncol = 2)
print(p3)

# ============================================================================
# 11. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Split UMAP by cell cycle phase
p4 <- DimPlot(pbmc, reduction = "umap", split.by = "Phase", ncol = 3)
print(p4)

# Violin plot showing cell cycle scores by cluster
p5 <- VlnPlot(pbmc, features = c("S.Score", "G2M.Score"), 
              group.by = "seurat_clusters", ncol = 2)
print(p5)

# ============================================================================
# 12. SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Number of clusters:", length(unique(pbmc$seurat_clusters)), "\n")
cat("\nCell cycle phase distribution:\n")
print(table(pbmc$Phase))
cat("\nCell cycle phase percentages:\n")
print(round(prop.table(table(pbmc$Phase)) * 100, 2))

# ============================================================================
# 13. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analyzed.rds")

# Save plots
ggsave("umap_cell_cycle.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("umap_clusters.png", plot = p2, width = 8, height = 6, dpi = 300)

cat("\nAnalysis complete! Results saved.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
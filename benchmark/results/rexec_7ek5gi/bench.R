setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# ============================================
# 1. LOAD DATA
# ============================================
# Assuming your data is in a standard format (e.g., 10X output)
# Adjust the path to your data directory
pbmc.data <- Read10X(data.dir = ".")

# Or if you have a count matrix:
# pbmc.data <- read.csv("your_count_matrix.csv", row.names = 1)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "PBMC_analysis", 
                           min.cells = 3,      # Keep genes expressed in >= 3 cells
                           min.features = 200) # Keep cells with >= 200 genes

# ============================================
# 2. QUALITY CONTROL
# ============================================
# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 2500 & 
                              percent.mt < 5)

# ============================================
# 3. NORMALIZATION
# ============================================
# Log-normalization (default method)
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Alternative: SCTransform (more advanced, recommended for newer analyses)
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# ============================================
# 4. IDENTIFY HIGHLY VARIABLE GENES
# ============================================
pbmc <- FindVariableFeatures(pbmc, 
                              selection.method = "vst", 
                              nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
print(paste("Top 10 variable genes:", paste(top10, collapse = ", ")))

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================
# 5. SCALING DATA
# ============================================
# Scale all genes
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Optional: regress out unwanted sources of variation
# pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt", "nCount_RNA"))

# ============================================
# 6. DIMENSIONALITY REDUCTION - PCA
# ============================================
pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# Heatmap of PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ============================================
# 7. DETERMINE SIGNIFICANT PCs
# ============================================
# JackStraw procedure (computationally intensive)
# pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
# JackStrawPlot(pbmc, dims = 1:20)

# ============================================
# 8. NON-LINEAR DIMENSIONALITY REDUCTION
# ============================================
# Choose number of PCs based on elbow plot (typically 10-30)
n_pcs <- 20

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

# Run t-SNE
pbmc <- RunTSNE(pbmc, dims = 1:n_pcs)

# Visualize
p1 <- DimPlot(pbmc, reduction = "umap") + ggtitle("UMAP")
p2 <- DimPlot(pbmc, reduction = "tsne") + ggtitle("t-SNE")
print(p1)
print(p2)

# ============================================
# 9. SAVE PROCESSED OBJECT
# ============================================
saveRDS(pbmc, file = "pbmc_processed.rds")

# Export normalized counts if needed
normalized_counts <- GetAssayData(pbmc, slot = "data")
write.csv(as.matrix(normalized_counts), "normalized_counts.csv")

# Export variable genes
variable_genes <- VariableFeatures(pbmc)
write.csv(variable_genes, "variable_genes.csv", row.names = FALSE)

# ============================================
# 10. SUMMARY STATISTICS
# ============================================
cat("\n=== Analysis Summary ===\n")
cat("Number of cells:", ncol(pbmc), "\n")
cat("Number of genes:", nrow(pbmc), "\n")
cat("Number of variable genes:", length(VariableFeatures(pbmc)), "\n")
cat("PCs used for UMAP/t-SNE:", n_pcs, "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
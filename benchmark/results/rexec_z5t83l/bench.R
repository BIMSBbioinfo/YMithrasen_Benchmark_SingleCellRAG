setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Single-cell RNA-seq Analysis Pipeline for 2,700 PBMCs
# Required package: Seurat

# Install Seurat if not already installed
if (!require("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

library(Seurat)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# 1. LOAD DATA
# ============================================================================
cat("Loading data...\n")

# Read the Matrix Market format data
pbmc_data <- ReadMtx(
  mtx = "matrix.mtx",
  cells = "barcodes.tsv",
  features = "genes.tsv",
  feature.column = 2  # Use HGNC symbols (column 2) as feature names
)

# Create Seurat object
pbmc <- CreateSeuratObject(
  counts = pbmc_data,
  project = "PBMC2700",
  min.cells = 3,      # Include features detected in at least 3 cells
  min.features = 200  # Include cells with at least 200 detected features
)

cat("Data loaded successfully!\n")
cat("Dimensions:", dim(pbmc), "\n")

# ============================================================================
# 2. QUALITY CONTROL (Optional but recommended)
# ============================================================================
cat("\nPerforming quality control...\n")

# Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
pdf("QC_plots.pdf", width = 12, height = 4)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", dim(pbmc), "\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================
cat("\nNormalizing data...\n")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ============================================================================
# 4. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================
cat("\nIdentifying highly variable features...\n")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable features:", paste(top10, collapse = ", "), "\n")

# Plot variable features
pdf("variable_features.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)
dev.off()

# ============================================================================
# 5. SCALE DATA
# ============================================================================
cat("\nScaling data...\n")

all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)

# ============================================================================
# 6. RUN PCA
# ============================================================================
cat("\nRunning PCA...\n")

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = 50
)

# Visualize PCA results
pdf("PCA_plots.pdf", width = 12, height = 8)
print(DimPlot(pbmc, reduction = "pca"))
print(ElbowPlot(pbmc, ndims = 50))
dev.off()

# Print top genes associated with PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# ============================================================================
# 7. CLUSTERING
# ============================================================================
cat("\nPerforming clustering...\n")

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# 8. RUN UMAP
# ============================================================================
cat("\nRunning UMAP...\n")

pbmc <- RunUMAP(pbmc, dims = 1:10)

# ============================================================================
# 9. VISUALIZATION
# ============================================================================
cat("\nGenerating visualizations...\n")

# UMAP plot colored by cluster
pdf("UMAP_clusters.pdf", width = 10, height = 8)
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6)
print(p1)
dev.off()

# UMAP plot with cluster labels and legend
pdf("UMAP_clusters_detailed.pdf", width = 12, height = 8)
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC 2,700 Cells - UMAP Clustering")
print(p2)
dev.off()

# Feature plots for marker genes
pdf("feature_plots.pdf", width = 12, height = 10)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================
cat("\nSaving results...\n")

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Export cluster assignments
cluster_assignments <- data.frame(
  barcode = colnames(pbmc),
  cluster = Idents(pbmc)
)
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

# Export UMAP coordinates
umap_coords <- as.data.frame(Embeddings(pbmc, reduction = "umap"))
umap_coords$barcode <- rownames(umap_coords)
umap_coords$cluster <- Idents(pbmc)
write.csv(umap_coords, "umap_coordinates.csv", row.names = FALSE)

# ============================================================================
# 11. SUMMARY STATISTICS
# ============================================================================
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(pbmc)))

cat("\nAnalysis complete! Check the generated PDF files for visualizations.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
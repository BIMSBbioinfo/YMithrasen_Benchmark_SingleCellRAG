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

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Read the Matrix Market format data
# Assumes files are in current working directory
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

cat("Loaded data dimensions:", dim(pbmc), "\n")

# ============================================================================
# 2. QUALITY CONTROL (Optional but recommended)
# ============================================================================

# Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", dim(pbmc), "\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

# Log-normalize the data
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("✓ Data normalized\n")

# ============================================================================
# 4. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================

# Find top 2000 highly variable genes
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable features:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

# ============================================================================
# 5. SCALE THE DATA
# ============================================================================

# Scale all genes
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)

cat("✓ Data scaled\n")

# ============================================================================
# 6. RUN PCA
# ============================================================================

# Perform linear dimensional reduction
pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc),
  npcs = 50
)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

cat("✓ PCA completed\n")

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Determine the dimensionality (using first 20 PCs based on typical PBMC data)
dims_use <- 1:20

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = dims_use)

# Find clusters (resolution 0.5 is good for ~3K cells)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check number of clusters
cat("Number of clusters identified:", length(levels(pbmc@active.ident)), "\n")
cat("Cluster sizes:\n")
print(table(Idents(pbmc)))

# ============================================================================
# 8. RUN UMAP
# ============================================================================

# Run UMAP for visualization
pbmc <- RunUMAP(pbmc, dims = dims_use)

cat("✓ UMAP completed\n")

# ============================================================================
# 9. VISUALIZE CLUSTERS
# ============================================================================

# UMAP plot colored by cluster
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) +
  ggtitle("PBMC Clusters (UMAP)")
print(p1)

# UMAP plot with cluster labels
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  NoLegend() +
  ggtitle("PBMC Clusters with Labels")
print(p2)

# PCA plot colored by cluster
p3 <- DimPlot(pbmc, reduction = "pca", label = TRUE) +
  ggtitle("PBMC Clusters (PCA)")
print(p3)

# ============================================================================
# 10. FIND CLUSTER MARKERS (Optional)
# ============================================================================

# Find markers for all clusters
pbmc_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Display top markers per cluster
top5_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top5_markers)

# Visualize top marker genes
top_genes <- top5_markers$gene[1:9]  # Top 9 genes for visualization
FeaturePlot(pbmc, features = top_genes, ncol = 3)

# Heatmap of top markers
top10_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Save cluster assignments
write.csv(
  data.frame(
    barcode = colnames(pbmc),
    cluster = Idents(pbmc)
  ),
  file = "cluster_assignments.csv",
  row.names = FALSE
)

# Save marker genes
write.csv(pbmc_markers, file = "cluster_markers.csv", row.names = FALSE)

cat("\n✓ Analysis complete! Results saved.\n")
cat("Seurat object: pbmc_analysis.rds\n")
cat("Cluster assignments: cluster_assignments.csv\n")
cat("Marker genes: cluster_markers.csv\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(levels(Idents(pbmc))), "\n")
cat("PCs used for clustering:", paste(dims_use, collapse = ", "), "\n")
cat("Clustering resolution:", 0.5, "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
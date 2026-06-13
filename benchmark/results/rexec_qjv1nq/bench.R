setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# ============================================================================
# PBMC scRNA-seq Analysis Pipeline with Seurat
# ============================================================================

# Install required packages (run once)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages <- c("Seurat", "dplyr", "ggplot2", "patchwork", "harmony")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
    }
}

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================================
# 1. DATA LOADING
# ============================================================================

# Option A: Load from 10X output directory
pbmc.data <- Read10X(data.dir = ".")

# Option B: Load from individual count matrix
# pbmc.data <- read.table("counts.txt", header = TRUE, row.names = 1)

# Create Seurat object
pbmc <- CreateSeuratObject(
    counts = pbmc.data,
    project = "PBMC_Analysis",
    min.cells = 3,      # Keep genes expressed in >= 3 cells
    min.features = 200  # Keep cells with >= 200 genes
)

# ============================================================================
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Calculate ribosomal percentage (optional)
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# Visualize QC metrics
VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.1)

# Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on QC metrics
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 5000 & 
                        percent.mt < 15)

cat("Cells after filtering:", ncol(pbmc), "\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

# Log-normalization
pbmc <- NormalizeData(pbmc, 
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

# Alternative: SCTransform (recommended for newer analyses)
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, 
                              selection.method = "vst",
                              nfeatures = 2000)

# Identify top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)
print(top10)

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# ============================================================================
# 5. SCALING AND PCA
# ============================================================================

# Scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, 
                  features = all.genes,
                  vars.to.regress = "percent.mt")  # Regress out MT%

# Run PCA
pbmc <- RunPCA(pbmc, 
               features = VariableFeatures(object = pbmc),
               npcs = 50)

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ============================================================================
# 6. DETERMINE DIMENSIONALITY
# ============================================================================

# Elbow plot
ElbowPlot(pbmc, ndims = 50)

# JackStraw (more rigorous but slower)
# pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
# JackStrawPlot(pbmc, dims = 1:50)

# Choose number of PCs (adjust based on elbow plot)
n_pcs <- 30

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)

# Find clusters (try different resolutions)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check cluster assignments
head(Idents(pbmc), 5)
table(Idents(pbmc))

# ============================================================================
# 8. NON-LINEAR DIMENSIONAL REDUCTION
# ============================================================================

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# Run tSNE (alternative)
pbmc <- RunTSNE(pbmc, dims = 1:n_pcs)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5)

# ============================================================================
# 9. FIND CLUSTER BIOMARKERS
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

# View top markers per cluster
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# Find markers for specific cluster (e.g., cluster 0)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)

# ============================================================================
# 10. VISUALIZATION OF MARKER GENES
# ============================================================================

# Violin plots
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# Feature plots on UMAP
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14",
                               "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# Heatmap of top markers
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()

# Dot plot
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP",
                     "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A",
                     "CD14", "LYZ", "S100A9",
                     "FCGR3A", "MS4A7",
                     "FCER1A", "CST3",
                     "PPBP")

DotPlot(pbmc, features = markers.to.plot, cols = c("blue", "red"), 
        dot.scale = 8) + RotatedAxis()

# ============================================================================
# 11. CELL TYPE ANNOTATION
# ============================================================================

# Manual annotation based on canonical markers
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
                     "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# ============================================================================
# 12. DIFFERENTIAL EXPRESSION BETWEEN CONDITIONS (if applicable)
# ============================================================================

# If you have multiple conditions/samples
# pbmc$condition <- factor(pbmc$orig.ident, levels = c("Control", "Treatment"))

# DE between conditions within a cell type
# Idents(pbmc) <- "condition"
# condition.markers <- FindMarkers(pbmc, ident.1 = "Treatment", ident.2 = "Control")

# ============================================================================
# 13. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Export marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)

# Export UMAP coordinates with metadata
umap_coords <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
umap_coords$cluster <- Idents(pbmc)
umap_coords$cell_type <- pbmc@active.ident
write.csv(umap_coords, "umap_coordinates.csv")

# Export normalized counts
normalized_counts <- GetAssayData(pbmc, slot = "data")
write.csv(as.matrix(normalized_counts), "normalized_counts.csv")

# ============================================================================
# 14. SESSION INFO
# ============================================================================

sessionInfo()

cat("\n=== Analysis Complete ===\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
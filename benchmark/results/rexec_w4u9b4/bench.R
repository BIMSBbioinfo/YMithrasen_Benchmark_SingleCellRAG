setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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

# Find markers for specific cluster
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)

# Visualize top markers
top_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top_markers$gene) + NoLegend()

# ============================================================================
# 10. CELL TYPE ANNOTATION
# ============================================================================

# Known PBMC markers
pbmc_markers <- c(
    "CD3D", "CD3E",           # T cells
    "CD4",                     # CD4+ T cells
    "CD8A", "CD8B",           # CD8+ T cells
    "IL7R", "CCR7",           # Naive T cells
    "GZMK",                    # Memory T cells
    "GNLY", "NKG7",           # NK cells
    "MS4A1", "CD79A",         # B cells
    "CD14", "LYZ",            # Monocytes
    "FCGR3A", "MS4A7",        # FCGR3A+ Monocytes
    "FCER1A", "CST3",         # Dendritic cells
    "PPBP"                     # Platelets
)

# Feature plots
FeaturePlot(pbmc, features = pbmc_markers[1:9], ncol = 3)

# Violin plots
VlnPlot(pbmc, features = pbmc_markers[1:6], ncol = 3)

# Dot plot
DotPlot(pbmc, features = pbmc_markers) + RotatedAxis()

# Assign cell type identities
new.cluster.ids <- c(
    "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", 
    "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet"
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# ============================================================================
# 11. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

# Compare two cell types
cd4_vs_cd8 <- FindMarkers(pbmc, 
                          ident.1 = "CD4 T", 
                          ident.2 = "CD8 T",
                          min.pct = 0.25)
head(cd4_vs_cd8, n = 10)

# ============================================================================
# 12. SAVE RESULTS
# ============================================================================

# Save Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Export marker genes
write.csv(pbmc.markers, "cluster_markers.csv", row.names = FALSE)

# Export metadata
write.csv(pbmc@meta.data, "cell_metadata.csv")

# Export UMAP coordinates
umap_coords <- Embeddings(pbmc, reduction = "umap")
write.csv(umap_coords, "umap_coordinates.csv")

# ============================================================================
# 13. ADDITIONAL VISUALIZATIONS
# ============================================================================

# Cluster composition
cluster_counts <- table(Idents(pbmc))
barplot(cluster_counts, las = 2, main = "Cells per Cluster")

# Feature plot for specific genes
FeaturePlot(pbmc, 
            features = c("CD3D", "CD14", "MS4A1", "NKG7"),
            ncol = 2)

# Ridge plot
RidgePlot(pbmc, features = c("CD3D", "CD14", "MS4A1"), ncol = 1)

# Generate summary report
cat("\n=== Analysis Summary ===\n")
cat("Total cells:", ncol(pbmc), "\n")
cat("Total genes:", nrow(pbmc), "\n")
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
cat("PCs used:", n_pcs, "\n")

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
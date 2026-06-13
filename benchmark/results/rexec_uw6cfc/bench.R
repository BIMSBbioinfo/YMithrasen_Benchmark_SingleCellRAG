setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required library
library(Seurat)

# Set seed for reproducibility
set.seed(42)

# 1. Load 10X Genomics data from working directory
pbmc.data <- Read10X(data.dir = ".")

# 2. Create Seurat object with quality filters
# Keep genes detected in at least 3 cells
# Keep cells with at least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          project = "PBMC_Analysis", 
                          min.cells = 3, 
                          min.features = 200)

# 3. Calculate mitochondrial QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4. Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 5. Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 6. Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 7. Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Plot variable features
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 8. Scale data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 9. Perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# 10. Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# 11. Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 12. Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE)

# 13. Find marker genes for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25,
                               test.use = "wilcox")

# Display top markers per cluster
top5_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
print(top5_markers)

# 14. Visualize top marker genes
VlnPlot(pbmc, features = head(top5_markers$gene, 9))
FeaturePlot(pbmc, features = head(top5_markers$gene, 9))

# 15. Generate heatmap of top markers
top10_markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10_markers$gene) + NoLegend()

# Save Seurat object
saveRDS(pbmc, file = "pbmc_analysis.rds")

# Export results
write.csv(pbmc.markers, file = "cluster_markers.csv", row.names = FALSE)
write.csv(pbmc@meta.data, file = "cell_metadata.csv")

print("Analysis complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)

# Read 10X data
pbmc.data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Find all markers
all.markers <- FindAllMarkers(pbmc, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Filter significant markers
significant.markers <- all.markers %>% filter(p_val_adj < 0.05)

# Get top 10 markers per cluster
top10.markers <- significant.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

print(top10.markers)

# Annotate cell types based on known markers
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Platelets")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# Visualize key marker genes
FeaturePlot(pbmc, features = c("CD3D", "CD3E", "CD79A", "MS4A1", "CD14", "LYZ", "NKG7", "GNLY"))

# Violin plot for key markers
VlnPlot(pbmc, features = c("CD3D", "CD79A", "CD14", "NKG7"))

# Heatmap of top markers
DoHeatmap(pbmc, features = top10.markers$gene) + NoLegend()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
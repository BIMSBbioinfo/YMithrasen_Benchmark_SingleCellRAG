setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)

# Read 10X data
pbmc_data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "PBMC", min.cells = 3, min.features = 200)

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

# Find marker genes for all clusters
pbmc_markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, test.use = "wilcox")

# Display top markers per cluster
top_markers <- pbmc_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
print(top_markers)

# Prepare data for SingleR annotation
pbmc_counts <- GetAssayData(pbmc, slot = "data")
pbmc_clusters <- pbmc$seurat_clusters

# Load reference dataset
ref <- MonacoImmuneData()

# Run SingleR for cell type annotation
singler_results <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main, clusters = pbmc_clusters)

# Add annotations to Seurat object
pbmc$celltype <- singler_results$labels[pbmc_clusters]

# Visualize clusters
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

# Visualize annotated cell types
DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)

# Display annotation results
print(table(pbmc$celltype, pbmc$seurat_clusters))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
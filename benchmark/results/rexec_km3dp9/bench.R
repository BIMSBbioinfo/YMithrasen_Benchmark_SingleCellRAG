setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required library
library(Seurat)

# Read 10x data from current directory
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
# Filter: keep genes expressed in ≥3 cells, keep cells with ≥200 genes
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          project = "PBMC", 
                          min.cells = 3, 
                          min.features = 200)

# Normalize data (log-normalization)
pbmc <- NormalizeData(pbmc, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
ElbowPlot(pbmc, ndims = 50)

# Find neighbors and clusters
# Using first 30 PCs (adjust based on ElbowPlot)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP with same dimensions
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Visualize clusters on UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# Optional: Save the plot
ggsave("pbmc_umap_clusters.png", width = 10, height = 8, dpi = 300)

# View cluster assignments
table(Idents(pbmc))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
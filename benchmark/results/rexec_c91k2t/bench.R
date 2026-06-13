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

# Load 10X Genomics data from current directory
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object with quality filtering
# min.cells = 3: Keep genes detected in at least 3 cells
# min.features = 200: Keep cells with at least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                          min.cells = 3, 
                          min.features = 200)

# Normalize the data
pbmc <- NormalizeData(pbmc)

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc)

# Determine dimensionality (examine elbow plot)
ElbowPlot(pbmc)

# Cluster the cells
# First compute the neighbor graph using first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP for visualization using first 10 PCs
pbmc <- RunUMAP(pbmc, dims = 1:10, reduction = "pca")

# Visualize clusters
DimPlot(pbmc, reduction = "umap")

# Print summary
print(pbmc)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
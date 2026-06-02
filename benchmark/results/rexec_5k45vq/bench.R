setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load necessary libraries
library(Seurat)
library(Matrix)
library(ggplot2)

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, rownames = genes[, 2], colnames = barcodes[, 1])

# Filter cells and genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 1000)
pbmc <- subset(pbmc, subset = nFeature_RNA > 0)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")

# Scale data
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine significant PCs
elbow_point <- WhichDimShouldIUse(pbmc, dims = 1:15)
pbmc <- RunUMAP(pbmc, dims = 1:elbow_point)

# Find clusters
pbmc <- FindNeighbors(pbmc, dims = 1:elbow_point)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Define marker genes for each cell type
marker_genes <- list(
  CD4_T = c("IL7R", "CCR7"),
  CD14_Monocytes = c("CD14", "LYZ"),
  B_cells = c("MS4A1"),
  CD8_T = c("CD8A"),
  NK_cells = c("GNLY", "NKG7")
)

# Annotate clusters with cell type labels
cluster_labels <- sapply(1:ncol(pbmc), function(x) {
  cluster <- pbmc@meta.data$seurat_clusters[x]
  max_score <- 0
  cell_type <- NA
  for (cell_type_name in names(marker_genes)) {
    score <- rowSums(pbmc@assays$RNA@data[which(genes[, 2] %in% marker_genes[[cell_type_name]]), x])
    if (score > max_score) {
      max_score <- score
      cell_type <- cell_type_name
    }
  }
  return(cell_type)
})

# Add cell type labels to Seurat object
pbmc@meta.data$cell_type <- cluster_labels

# Visualize UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", group.by = "cell_type", cols = "coolwarm")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
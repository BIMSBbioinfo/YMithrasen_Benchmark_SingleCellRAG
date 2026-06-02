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

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, rownames = genes[, 2], colnames = barcodes[, 1])

# Filter cells and genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst")

# Scale data
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine significant PCs
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
elbow <- which.max(pbmc@reductions$pca@jackstraw$p.values[, 2])
pbmc <- RunUMAP(pbmc, dims = 1:elbow)

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:elbow)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Define marker genes for each cell type
marker_genes <- list(
  "CD4 T cells" = c("IL7R", "CCR7"),
  "CD14 Monocytes" = c("CD14", "LYZ"),
  "B cells" = c("MS4A1"),
  "CD8 T cells" = c("CD8A"),
  "NK cells" = c("GNLY", "NKG7")
)

# Find markers for each cluster
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Annotate clusters with cell type labels
cell_type_labels <- character(nrow(pbmc@meta.data))
for (i in 1:ncol(cluster_markers)) {
  cluster <- cluster_markers[[i]]
  for (cell_type in names(marker_genes)) {
    markers <- marker_genes[[cell_type]]
    if (any(cluster$gene %in% markers)) {
      cell_type_labels[pbmc@meta.data$seurat_clusters == i] <- cell_type
    }
  }
}

# Add cell type labels to Seurat object
pbmc@meta.data$cell_type <- cell_type_labels

# Visualize UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", group.by = "cell_type")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(Matrix)
library(ggplot2)

# Load data
mat <- readMM("matrix.mtx")
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)$V1
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
rownames(mat) <- genes$V2   # use HGNC symbols
colnames(mat) <- barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = mat, project = "PBMC", min.cells = 3, min.features = 200)

# QC: calculate percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization and variable features
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling and PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Identify markers for each cluster
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Assign cell type labels based on known markers
cluster_ids <- c(
  "CD4 T cells" = which(markers$gene %in% c("IL7R", "CCR7") & markers$p_val_adj < 0.05)[1],
  "CD14 Monocytes" = which(markers$gene %in% c("CD14", "LYZ") & markers$p_val_adj < 0.05)[1],
  "B cells" = which(markers$gene %in% c("MS4A1") & markers$p_val_adj < 0.05)[1],
  "CD8 T cells" = which(markers$gene %in% c("CD8A") & markers$p_val_adj < 0.05)[1],
  "NK cells" = which(markers$gene %in% c("GNLY", "NKG7") & markers$p_val_adj < 0.05)[1]
)

# Map cluster numbers to labels
cluster_levels <- levels(pbmc)
cluster_labels <- setNames(names(cluster_ids), cluster_ids)
pbmc <- RenameIdents(pbmc, cluster_labels)

# Plot annotated UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Annotated PBMC UMAP")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
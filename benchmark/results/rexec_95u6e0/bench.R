setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
counts <- readMM("matrix.mtx")
rownames(counts) <- genes$V2   # use HGNC symbols
colnames(counts) <- barcodes$V1

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC")

# QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization and feature selection
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling, PCA, clustering
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Identify markers
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Manual annotation based on marker genes
cluster_mapping <- c(
  "0" = "CD4 T cells",
  "1" = "CD14 Monocytes",
  "2" = "B cells",
  "3" = "CD8 T cells",
  "4" = "NK cells"
)
pbmc$celltype <- factor(cluster_mapping[as.character(pbmc$seurat_clusters)], 
                        levels = c("CD4 T cells", "CD14 Monocytes", "B cells", "CD8 T cells", "NK cells"))

# Plot annotated UMAP
DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  ggtitle("Annotated PBMC UMAP")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(Matrix)
library(dplyr)

# Load data
matrix_file <- "matrix.mtx"
barcodes_file <- "barcodes.tsv"
genes_file <- "genes.tsv"

counts <- readMM(matrix_file)
barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- genes$V2  # HGNC symbols
colnames(counts) <- barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Standard workflow
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Find markers
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get top 10 markers per cluster
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Heatmap of top markers
DoHeatmap(pbmc, features = top10$gene, group.by = "cluster", nfeatures = 10, 
          size = 3, angle = 45) + 
  theme_minimal()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
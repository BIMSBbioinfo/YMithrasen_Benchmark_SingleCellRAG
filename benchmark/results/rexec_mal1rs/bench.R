setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)
library(DoubletFinder)
library(ggplot2)

# Set seed for reproducibility
set.seed(1234)

# Read data
matrix_file <- "matrix.mtx"
barcodes_file <- "barcodes.tsv"
genes_file <- "genes.tsv"

counts <- readMM(matrix_file)
barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1
genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- genes$V1
colnames(counts) <- barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Basic preprocessing
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# DoubletFinder parameters
pbmc <- doubletFinder_v3(pbmc, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = 0.1 * ncol(pbmc), reuse.pANN = FALSE, pANN = 0.25, PCs = 1:30, pK = 0.09, nExp = 0.1 * ncol(pbmc), reuse.pANN = FALSE, pANN = 0.25)

# Add doublet classification to metadata
pbmc$doublet_class <- ifelse(pbmc$DF.classifications_0.25_0.09_0.1 == "Doublet", "Doublet", "Singlet")

# Visualize UMAP colored by doublet status
DimPlot(pbmc, reduction = "umap", group.by = "doublet_class", pt.size = 0.5) + 
  ggtitle("UMAP of PBMCs: Singlets vs Doublets")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
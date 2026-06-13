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

# QC metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine dimensionality
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = 0.1 * ncol(pbmc), reuse.pANN = FALSE, sct = FALSE)

# Visualize singlets vs doublets
pbmc$doublet_class <- ifelse(pbmc$DF.classifications_0.25_0.09_0.1 == "Doublet", "Doublet", "Singlet")
DimPlot(pbmc, reduction = "umap", group.by = "doublet_class", pt.size = 0.5) + ggtitle("UMAP: Singlets vs Doublets")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
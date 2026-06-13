setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(DoubletFinder)
  library(ggplot2)
})

# Read data
barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes   <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
counts  <- readMM("matrix.mtx")

# Assign row and column names
rownames(counts) <- genes$V1          # gene IDs
colnames(counts) <- barcodes$V1       # cell barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add gene symbols as metadata
pbmc[["gene_symbols"]] <- genes$V2[match(rownames(pbmc), genes$V1)]

# Standard Seurat preprocessing
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# DoubletFinder parameters
pbmc <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE, nExp = 0.25, PCs = 1:30)
sweep_res <- summarizeSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
bcmvn <- find.pK(pbmc, PCs = 1:30, sct = FALSE)

# Choose pK (here we pick the median of the top 5 pK values)
pK <- median(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:30, pN = 0.25, pK = pK, nExp = 0.25, reuse.pANN = FALSE, sct = FALSE)

# Add doublet classification to metadata
pbmc$doublet_class <- pbmc$DF.classifications_0.25_0.25_0.25

# Visualize UMAP colored by doublet status
DimPlot(pbmc, reduction = "umap", group.by = "doublet_class", pt.size = 0.5) + 
  ggtitle("UMAP: Singlets vs Doublets")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
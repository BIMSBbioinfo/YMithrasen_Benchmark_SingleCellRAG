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

# Create a Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC", 
                           min.cells = 3, 
                           min.features = 3)

# Add cell and gene metadata
pbmc@colData <- data.frame(barcodes = barcodes$V1)
pbmc@rowData <- data.frame(genes = genes$V1, 
                           symbols = genes$V2)

# Calculate the number of detected genes per cell
pbmc <- CalculateSizeFactors(object = pbmc)

# Calculate the percentage of mitochondrial genes
mito.genes <- genes[grep("^MT-", genes$V2), ]
pbmc <- AddModuleScore(object = pbmc, 
                        features = mito.genes$V1, 
                        name = "percent.mt")

# Calculate the total UMI counts per cell
pbmc <- CalculateSizeFactors(object = pbmc)

# Filter out low-quality cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize the data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
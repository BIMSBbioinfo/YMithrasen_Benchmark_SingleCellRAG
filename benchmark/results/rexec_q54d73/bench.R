setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load necessary libraries
library(Seurat)
library(ggplot2)

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, 
                                    project = "PBMCs", 
                                    min.cells = 3, 
                                    min.features = 200)

# Add cell and gene metadata
seurat_object@colData <- data.frame(barcodes = barcodes$V1)
seurat_object@rowData <- data.frame(genes = genes$V1, 
                                      symbols = genes$V2)

# Normalize data
seurat_object <- NormalizeData(seurat_object, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# Identify highly variable features
seurat_object <- FindVariableFeatures(seurat_object, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)

# Scale data
seurat_object <- ScaleData(seurat_object, 
                            features = VariableFeatures(object = seurat_object))

# Run PCA
seurat_object <- RunPCA(seurat_object, 
                         features = VariableFeatures(object = seurat_object))

# Perform clustering
seurat_object <- FindNeighbors(seurat_object, 
                               dims = 1:10)
seurat_object <- FindClusters(seurat_object, 
                              resolution = 0.5)

# Compute UMAP
seurat_object <- RunUMAP(seurat_object, 
                         dims = 1:10)

# Visualize clusters
DimPlot(seurat_object, 
         reduction = "umap", 
         group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
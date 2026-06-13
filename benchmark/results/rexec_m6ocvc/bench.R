setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(DoubletFinder)
library(scDblFinder)

# Set working directory and file paths
setwd("path/to/your/data")
sample_files <- c("sample1_filtered_feature_bc_matrix",
                  "sample2_filtered_feature_bc_matrix",
                  "sample3_filtered_feature_bc_matrix")

# Read 10x data and merge into one Seurat object
data_list <- lapply(sample_files, Read10X)
seurat_obj <- merge(data_list[[1]], y = data_list[-1], add.cell.ids = c("S1","S2","S3"))

# Basic QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Filter cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization and variable feature selection
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling and PCA
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Determine dimensionality
ElbowPlot(seurat_obj)

# UMAP and clustering
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Doublet detection with DoubletFinder
seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = 0.25 * ncol(seurat_obj), reuse.pANN = FALSE, sct = FALSE)
# Identify doublet cluster(s)
doublet_clusters <- WhichCells(seurat_obj, idents = "DF_class_1")  # adjust if multiple classes
seurat_obj <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), doublet_clusters))

# Alternatively, doublet detection with scDblFinder
seurat_obj <- scDblFinder(seurat_obj, pN = 0.25, pK = 0.09, nExp = 0.25 * ncol(seurat_obj), sct = FALSE)
# Remove predicted doublets
seurat_obj <- subset(seurat_obj, subset = scDblFinder.classification != "Doublet")

# Proceed with downstream analysis
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
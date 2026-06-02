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
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC", 
                           min.cells = 3, 
                           min.features = 200)

# Add gene and cell metadata
pbmc@colnames <- barcodes$V1
pbmc@rownames <- genes$V1

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

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Run clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Annotate clusters with known cell type labels
cd4_t_cells <- WhichCells(pbmc, expression = rowMeans(pbmc@assays$RNA@data[rownames(pbmc@assays$RNA@data) %in% c("IL7R", "CCR7"), ]))
cd14_monocytes <- WhichCells(pbmc, expression = rowMeans(pbmc@assays$RNA@data[rownames(pbmc@assays$RNA@data) %in% c("CD14", "LYZ"), ]))
b_cells <- WhichCells(pbmc, expression = rowMeans(pbmc@assays$RNA@data[rownames(pbmc@assays$RNA@data) %in% c("MS4A1"), ]))
cd8_t_cells <- WhichCells(pbmc, expression = rowMeans(pbmc@assays$RNA@data[rownames(pbmc@assays$RNA@data) %in% c("CD8A"), ]))
nk_cells <- WhichCells(pbmc, expression = rowMeans(pbmc@assays$RNA@data[rownames(pbmc@assays$RNA@data) %in% c("GNLY", "NKG7"), ]))

cell_types <- rep("Unknown", ncol(pbmc))
cell_types[cd4_t_cells] <- "CD4 T cells"
cell_types[cd14_monocytes] <- "CD14 Monocytes"
cell_types[b_cells] <- "B cells"
cell_types[cd8_t_cells] <- "CD8 T cells"
cell_types[nk_cells] <- "NK cells"

# Update Seurat object with cell type labels
pbmc@meta.data$cell_type <- cell_types

# Visualize UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", group.by = "cell_type")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
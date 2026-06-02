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
                           min.features = 200,
                           cell.names = barcodes$V1,
                           gene.names = genes$V1)

# Filter cells and genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA > 1000)

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
elbow <- WhichDimAreSignificant(pbmc, dims = 1:20)
pbmc <- RunUMAP(pbmc, dims = 1:elbow)

# Find clusters
pbmc <- FindNeighbors(pbmc, dims = 1:elbow)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Define marker genes for each cell type
marker_genes <- list(
  CD4_T = c("IL7R", "CCR7"),
  CD14_Monocytes = c("CD14", "LYZ"),
  B_cells = c("MS4A1"),
  CD8_T = c("CD8A"),
  NK_cells = c("GNLY", "NKG7")
)

# Find markers for each cluster
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Annotate clusters with cell type labels
cell_type_labels <- sapply(1:pbmc@meta.data$seurat_clusters, function(x) {
  cluster_marker <- cluster_markers %>% filter(cluster == x)
  max_score <- 0
  cell_type <- NA
  for (cell_type in names(marker_genes)) {
    score <- sum(cluster_marker %>% filter(gene %in% marker_genes[[cell_type]]) %>% pull(p_val_adj))
    if (score > max_score) {
      max_score <- score
      cell_type <- cell_type
    }
  }
  return(cell_type)
})

# Add cell type labels to Seurat object
pbmc@meta.data$cell_type <- cell_type_labels[pbmc@meta.data$seurat_clusters]

# Visualize UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", group.by = "cell_type")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
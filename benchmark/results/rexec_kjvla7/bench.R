setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)

# Read data
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
counts <- readMM("matrix.mtx")

# Assign row and column names
rownames(counts) <- genes$V1          # gene IDs
colnames(counts) <- barcodes$V1       # cell barcodes

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add gene symbols as metadata
pbmc[["gene_symbols"]] <- genes$V2

# QC: filter cells and genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization and identification of variable features
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling and PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine dimensionality
ElbowPlot(pbmc)

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Identify markers for each cluster
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Function to assign cell type based on marker genes
assign_celltype <- function(cluster_id, markers_df) {
  cluster_markers <- markers_df[markers_df$cluster == cluster_id, ]
  genes_present <- cluster_markers$gene
  
  if (any(c("IL7R", "CCR7") %in% genes_present)) {
    return("CD4 T cells")
  } else if (any(c("CD14", "LYZ") %in% genes_present)) {
    return("CD14 Monocytes")
  } else if (any("MS4A1" %in% genes_present)) {
    return("B cells")
  } else if (any("CD8A" %in% genes_present)) {
    return("CD8 T cells")
  } else if (any(c("GNLY", "NKG7") %in% genes_present)) {
    return("NK cells")
  } else {
    return("Other")
  }
}

# Apply assignment to all clusters
cluster_ids <- levels(pbmc$seurat_clusters)
celltype_labels <- sapply(cluster_ids, assign_celltype, markers_df = markers)
names(celltype_labels) <- cluster_ids

# Add cell type metadata
pbmc$celltype <- factor(celltype_labels[pbmc$seurat_clusters])

# Plot annotated UMAP
DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 4) + 
  ggtitle("Annotated UMAP of PBMCs") + 
  theme_minimal()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
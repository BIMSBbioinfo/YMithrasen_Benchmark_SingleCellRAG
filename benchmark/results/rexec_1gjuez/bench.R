setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)

# Load 10X data from current directory
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC", min.cells = 3, min.features = 200)

# Standard preprocessing workflow
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Find marker genes for each cluster
Idents(pbmc) <- pbmc$seurat_clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

# Filter significant markers
pbmc.markers.sig <- pbmc.markers %>% filter(p_val_adj < 0.05)

# Get top markers per cluster
top_markers <- pbmc.markers.sig %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top_markers)

# Automated cell type annotation using SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
counts_matrix <- GetAssayData(pbmc, slot = "data")
singler_results <- SingleR(test = counts_matrix, ref = ref, labels = ref$label.main)

# Add annotations to Seurat object
pbmc$celltype <- singler_results$labels

# Visualize clusters and cell types
DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)

# Summary of cell type annotations
table(pbmc$celltype, pbmc$seurat_clusters)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
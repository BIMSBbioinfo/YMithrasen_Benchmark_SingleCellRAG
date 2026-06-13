setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load your Seurat object (replace "pbmc_seurat.rds" with your file)
pbmc <- readRDS("pbmc_seurat.rds")

# Ensure clusters are identified (if not, run clustering)
if (!"seurat_clusters" %in% colnames(pbmc@meta.data)) {
  pbmc <- FindNeighbors(pbmc, dims = 1:20)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
}

# Select two clusters for comparison (e.g., cluster 0 vs cluster 1)
cluster1 <- "0"
cluster2 <- "1"

# Subset cells belonging to the two clusters
cells_to_compare <- WhichCells(pbmc, idents = c(cluster1, cluster2))
pbmc_subset <- subset(pbmc, cells = cells_to_compare)

# Perform differential expression analysis
de_markers <- FindMarkers(pbmc_subset, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = 0.25)

# Add gene symbols as a column if not present
de_markers <- de_markers %>% rownames_to_column(var = "gene")

# Visualize top 10 markers with violin plots
top_genes <- de_markers %>% arrange(p_val_adj) %>% slice_head(n = 10) %>% pull(gene)

# Violin plots
VlnPlot(pbmc_subset, features = top_genes, group.by = "seurat_clusters", pt.size = 0.1) + 
  ggtitle(paste("Top 10 DE genes between cluster", cluster1, "and cluster", cluster2))

# Feature plots for the top 3 markers
FeaturePlot(pbmc_subset, features = top_genes[1:3], cols = c("lightgrey", "blue")) + 
  ggtitle(paste("Feature plots for top 3 DE genes between cluster", cluster1, "and cluster", cluster2))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
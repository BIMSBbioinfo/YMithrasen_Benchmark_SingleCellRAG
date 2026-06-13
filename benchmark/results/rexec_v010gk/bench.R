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

# Load your PBMC dataset (replace with your actual data loading step)
# Example: pbmc <- readRDS("pbmc_seurat_object.rds")
pbmc <- readRDS("pbmc_seurat_object.rds")

# Ensure clusters are identified (if not, run clustering)
if (!"seurat_clusters" %in% colnames(pbmc@meta.data)) {
  pbmc <- FindNeighbors(pbmc, dims = 1:20)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
}

# Specify the two clusters of interest
cluster1 <- "0"
cluster2 <- "1"

# Subset cells belonging to the two clusters
cells_to_compare <- WhichCells(pbmc, idents = c(cluster1, cluster2))
pbmc_subset <- subset(pbmc, cells = cells_to_compare)

# Perform differential expression analysis
de_results <- FindMarkers(pbmc_subset, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = 0.25)

# Add gene names as a column
de_results$gene <- rownames(de_results)

# Visualize top differentially expressed genes with a violin plot
top_genes <- de_results %>% arrange(p_val_adj) %>% slice_head(n = 10) %>% pull(gene)

# Violin plots for top genes
pbmc_subset %>%
  VlnPlot(features = top_genes, group.by = "seurat_clusters", pt.size = 0.1) +
  theme_bw()

# Alternatively, create a heatmap of top genes
top_genes_matrix <- FetchData(pbmc_subset, vars = top_genes)
heatmap_data <- t(scale(t(top_genes_matrix)))
pheatmap::pheatmap(heatmap_data,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   main = paste("Top DE genes between clusters", cluster1, "and", cluster2))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
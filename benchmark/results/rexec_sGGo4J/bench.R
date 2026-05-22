setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))

# Check current cluster identities
current_clusters <- levels(Idents(pbmc))
print(paste("Current clusters:", paste(current_clusters, collapse = ", ")))

# Create named vector mapping cluster numbers to cell types
# Adjust these mappings based on your marker gene expression analysis
new_cluster_ids <- c(
  "Naive CD4 T",
  "CD14+ Mono",
  "Memory CD4 T",
  "B",
  "CD8 T",
  "FCGR3A+ Mono",
  "NK",
  "DC",
  "Platelet"
)

# Ensure the vector length matches number of clusters
if(length(new_cluster_ids) != length(current_clusters)) {
  # Adjust based on actual number of clusters
  new_cluster_ids <- new_cluster_ids[1:length(current_clusters)]
}

names(new_cluster_ids) <- current_clusters

# Rename cluster identities
pbmc <- RenameIdents(pbmc, new_cluster_ids)

# Generate UMAP plot with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# Alternative with legend
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5, repel = TRUE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
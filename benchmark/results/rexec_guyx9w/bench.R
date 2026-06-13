setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)

# Load your Seurat object (replace "pbmc_seurat.rds" with your file)
pbmc <- readRDS("pbmc_seurat.rds")

# Find markers for each cluster
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers to a CSV
write.csv(markers, file = "pbmc_cluster_markers.csv", row.names = FALSE)

# Annotate clusters using known marker genes
# Define a simple annotation dictionary
cluster_annotations <- list(
  "0" = "CD14+ Monocytes",
  "1" = "CD4 T cells",
  "2" = "B cells",
  "3" = "NK cells",
  "4" = "FCGR3A+ Monocytes",
  "5" = "Dendritic cells",
  "6" = "Megakaryocytes",
  "7" = "CD8 T cells",
  "8" = "Proliferating cells",
  "9" = "Naive B cells"
)

# Add annotations to metadata
pbmc$cell_type <- factor(cluster_annotations[as.character(pbmc$seurat_clusters)], levels = cluster_annotations)

# Save the annotated Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
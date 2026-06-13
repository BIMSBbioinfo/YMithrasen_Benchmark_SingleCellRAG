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

# View top markers per cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Basic cell type annotation using known markers
celltype_markers <- list(
  "CD4 T" = c("IL7R", "CCR7"),
  "CD8 T" = c("NKG7", "GZMB"),
  "NK" = c("NCAM1", "KLRD1"),
  "B" = c("MS4A1", "CD79A"),
  "Monocyte" = c("LYZ", "S100A8"),
  "Dendritic" = c("FCER1A", "CLEC10A"),
  "Megakaryocyte" = c("PPBP", "PF4")
)

# Assign cell type based on marker expression
pbmc$cell_type <- NA
for (ct in names(celltype_markers)) {
  pbmc$cell_type[which(
    rowSums(pbmc@assays$RNA@data[celltype_markers[[ct]], ] > 0) > 0
  )] <- ct
}

# Add annotations to metadata
pbmc@meta.data$cell_type <- pbmc$cell_type

# Save updated Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
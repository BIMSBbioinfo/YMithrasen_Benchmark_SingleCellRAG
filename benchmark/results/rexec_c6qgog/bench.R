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
write.csv(markers, "pbmc_cluster_markers.csv", row.names = FALSE)

# Annotate clusters using known marker genes
# Define a simple annotation dictionary
annotation_dict <- list(
  "CD3D" = "T cells",
  "MS4A1" = "B cells",
  "GNLY" = "NK cells",
  "FCER1A" = "Dendritic cells",
  "PPBP" = "Platelets",
  "LYZ" = "Monocytes"
)

# Function to assign cell type based on top marker genes
assign_cell_type <- function(cluster_id, markers_df, dict) {
  top_genes <- markers_df %>% filter(cluster == cluster_id) %>% arrange(desc(avg_log2FC)) %>% head(10) %>% pull(gene)
  for (gene in top_genes) {
    if (gene %in% names(dict)) {
      return(dict[[gene]])
    }
  }
  return("Unknown")
}

# Apply annotation
cluster_ids <- unique(markers$cluster)
cluster_annotations <- sapply(cluster_ids, assign_cell_type, markers_df = markers, dict = annotation_dict)

# Add annotations to Seurat object metadata
pbmc$cell_type <- factor(cluster_annotations[as.character(pbmc$seurat_clusters)])

# Save updated Seurat object
saveRDS(pbmc, "pbmc_seurat_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
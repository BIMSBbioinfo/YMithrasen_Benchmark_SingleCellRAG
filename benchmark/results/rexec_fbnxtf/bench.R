setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)

# Set seed for reproducibility
set.seed(1234)

# Define file paths
count_file <- "counts_matrix.csv"          # Path to raw count matrix (genes x cells)
output_dir <- "Seurat_analysis"           # Main output directory
plots_dir  <- file.path(output_dir, "plots")
markers_dir <- file.path(output_dir, "markers")

# Create output directories if they don't exist
dir.create(output_dir, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)
dir.create(markers_dir, showWarnings = FALSE)

# Read raw count data
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC", min.cells = 3, min.features = 200)

# Add mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# SCTransform normalization
pbmc <- SCTransform(pbmc, verbose = FALSE)

# Run PCA
pbmc <- RunPCA(pbmc, verbose = FALSE)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

# Save Seurat object
saveRDS(pbmc, file = file.path(output_dir, "pbmc_seurat_object.rds"))
SaveH5Seurat(pbmc, filename = file.path(output_dir, "pbmc_seurat_object.h5Seurat"))

# UMAP plot
umap_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of PBMC")
ggsave(filename = file.path(plots_dir, "UMAP_plot.png"), plot = umap_plot, width = 8, height = 6, dpi = 300)

# Find cluster markers
cluster_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers to CSV
write.csv(cluster_markers, file = file.path(markers_dir, "cluster_markers.csv"), row.names = FALSE)

# Save top markers per cluster
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top_markers, file = file.path(markers_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)

# Save Seurat object as RDS for downstream use
saveRDS(pbmc, file = file.path(output_dir, "pbmc_seurat_final.rds"))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
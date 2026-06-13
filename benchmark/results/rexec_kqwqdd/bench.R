setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)

# Set seed for reproducibility
set.seed(42)

# Read 10X data
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# Run PCA
pbmc <- RunPCA(pbmc, npcs = 50)

# Elbow plot to determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)

# Find all markers
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

print(top_markers)

# Automated cell type annotation using SingleR
ref <- celldex::HumanPrimaryCellAtlasData()
pbmc_counts <- GetAssayData(pbmc, slot = "data")
pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)

# Add predictions to Seurat object
pbmc[["SingleR.labels"]] <- pred$labels

# Visualize annotated cell types
DimPlot(pbmc, reduction = "umap", group.by = "SingleR.labels", label = TRUE)

# Save results
saveRDS(pbmc, file = "pbmc_analyzed.rds")
write.csv(pbmc.markers, file = "pbmc_markers.csv", row.names = FALSE)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
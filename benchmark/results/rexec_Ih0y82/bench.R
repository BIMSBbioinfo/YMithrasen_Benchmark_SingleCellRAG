setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(patchwork)))

# Load 10X data from current directory
sc.data <- Read10X(data.dir = ".")

# Create Seurat object
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# QC filtering
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Scale data
sc <- ScaleData(sc, features = rownames(sc))

# Run PCA
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

# Find neighbors and clusters
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

# Run UMAP
sc <- RunUMAP(sc, dims = 1:10)

# Find all markers
all.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Filter by significance and select top 5 per cluster
top5 <- all.markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

# Create heatmap
DoHeatmap(sc, features = top5$gene) + NoLegend()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
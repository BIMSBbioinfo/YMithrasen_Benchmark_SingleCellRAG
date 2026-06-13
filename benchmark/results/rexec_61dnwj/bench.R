setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(ggplot2)

# Load data
sc.data <- Read10X(data.dir = ".")
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# QC metrics
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# QC filtering
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Scale data
sc <- ScaleData(sc, features = rownames(sc))

# PCA
sc <- RunPCA(sc, features = VariableFeatures(object = sc))

# Clustering
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

# UMAP
sc <- RunUMAP(sc, dims = 1:10)

# DoubletFinder parameter sweep
sweep.res.list <- paramSweep_v3(sc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Extract optimal pK
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Estimate doublet proportion
annotations <- sc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
sc <- doubletFinder_v3(sc, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Find DoubletFinder classification column
df_col <- grep("^DF.classifications", colnames(sc@meta.data), value = TRUE)

# Visualize doublets
DimPlot(sc, reduction = "umap", group.by = df_col)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
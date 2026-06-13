setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(DoubletFinder)
library(ggplot2)

# Load 10X data
sc.data <- Read10X(data.dir = ".")

# Create Seurat object
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Add mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# Preprocessing pipeline
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc, features = VariableFeatures(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc))
sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:10)

# DoubletFinder parameter sweep
sweep.res.list <- paramSweep_v3(sc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Use optimal pK or fixed value
pK <- 0.09

# Calculate expected doublets
annotations <- sc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
sc <- doubletFinder_v3(sc, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get classification column name
classification_col <- grep("^DF.classifications", colnames(sc@meta.data), value = TRUE)

# Visualize doublets on UMAP
DimPlot(sc, reduction = "umap", group.by = classification_col)

# Subset to singlets only
sc_singlet <- subset(sc, subset = !!sym(classification_col) == "Singlet")

# Visualize singlets
DimPlot(sc_singlet, reduction = "umap", group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
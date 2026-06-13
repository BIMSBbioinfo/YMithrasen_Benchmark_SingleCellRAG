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
library(Matrix)

# Load 10X data
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc2700", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# DoubletFinder - pK identification
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select optimal pK
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Homotypic doublet proportion estimate
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# Expected number of doublets (assuming 7% doublet rate)
nExp_poi <- round(0.07 * nrow(pbmc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get the doublet classification column name
doublet_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets on UMAP
DimPlot(pbmc, reduction = "umap", group.by = doublet_col, pt.size = 0.5) + 
  ggtitle("Doublet Detection Results")

# Subset singlets only
pbmc_singlet <- subset(pbmc, subset = !!sym(doublet_col) == "Singlet")

# Visualize singlets only
DimPlot(pbmc_singlet, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5) + 
  ggtitle("Singlets Only")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
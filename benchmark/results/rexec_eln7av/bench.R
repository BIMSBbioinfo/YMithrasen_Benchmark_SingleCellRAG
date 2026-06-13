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
pbmc_data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "PBMC", min.cells = 3, min.features = 200)

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Standard preprocessing
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# pK identification
sweep.res.list <- paramSweep(pbmc, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot pK optimization
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# Use optimal pK (or fixed value)
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected doublets
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.07 * ncol(pbmc))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder(pbmc, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get classification column name
df_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets
DimPlot(pbmc, reduction = "umap", group.by = df_col)

# Filter to keep only singlets
pbmc_singlet <- subset(pbmc, subset = pbmc@meta.data[[df_col]] == "Singlet")

# Visualize singlets
DimPlot(pbmc_singlet, reduction = "umap", group.by = "seurat_clusters")

# Print summary
table(pbmc@meta.data[[df_col]])
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
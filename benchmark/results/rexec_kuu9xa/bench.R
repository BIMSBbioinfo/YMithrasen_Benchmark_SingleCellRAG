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

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize and find variable features
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data and run PCA
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

# Run UMAP and find neighbors/clusters
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# pK identification using paramSweep
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot BCmetric to identify optimal pK
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# Use optimal pK (typically where BCmetric peaks)
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected doublets
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(pbmc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get the classification column name
df_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets
DimPlot(pbmc, reduction = "umap", group.by = df_col)

# Filter to keep only singlets
pbmc_singlets <- subset(pbmc, subset = pbmc@meta.data[[df_col]] == "Singlet")

# Visualize singlets only
DimPlot(pbmc_singlets, reduction = "umap", group.by = "seurat_clusters")

# Print summary
cat("Original cells:", ncol(pbmc), "\n")
cat("Singlets retained:", ncol(pbmc_singlets), "\n")
cat("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
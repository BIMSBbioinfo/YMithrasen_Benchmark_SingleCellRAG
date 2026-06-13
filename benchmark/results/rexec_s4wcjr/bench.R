setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(DoubletFinder)

# Load 10X data
pbmc.data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

# Standard preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# DoubletFinder parameter optimization
sweep.res <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected doublets (0.8% per 1000 cells)
nExp_poi <- round(0.008 * ncol(pbmc) * (ncol(pbmc)/1000))

# Estimate homotypic doublets
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get classification column name
df_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Filter to keep only singlets
pbmc_singlets <- subset(pbmc, subset = !!sym(df_col) == "Singlet")

# Summary
cat("Total cells:", ncol(pbmc), "\n")
cat("Doublets detected:", sum(pbmc@meta.data[[df_col]] == "Doublet"), "\n")
cat("Singlets retained:", ncol(pbmc_singlets), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
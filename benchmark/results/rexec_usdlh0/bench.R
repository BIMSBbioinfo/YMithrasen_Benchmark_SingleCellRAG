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

# Load 10X data from current directory
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

# Standard preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)

# Dimensionality reduction (required for DoubletFinder)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Find clusters for homotypic doublet estimation
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Optimize pK parameter using paramSweep
sweep.res <- paramSweep_v3(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select optimal pK
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected doublets (7.5% rate for ~10K cells)
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * ncol(pbmc))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = pK_optimal, 
                         nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get classification column name (it's dynamically generated)
classification_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets
DimPlot(pbmc, reduction = "umap", group.by = classification_col)

# Subset to keep only singlets
pbmc_singlets <- subset(pbmc, subset = pbmc@meta.data[[classification_col]] == "Singlet")

# Show summary
cat("Original cells:", ncol(pbmc), "\n")
cat("Singlets:", ncol(pbmc_singlets), "\n")
cat("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets), "\n")

# Visualize singlets only
DimPlot(pbmc_singlets, reduction = "umap", group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
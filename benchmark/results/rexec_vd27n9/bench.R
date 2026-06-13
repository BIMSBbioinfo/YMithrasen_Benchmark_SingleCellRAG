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

# Add mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Standard preprocessing workflow
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# Run PCA and UMAP (required for DoubletFinder)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:20)

# Find clusters for homotypic doublet estimation
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Optimize pK parameter using paramSweep
sweep.res <- paramSweep_v3(pbmc, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot pK optimization
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal()

# Select optimal pK (highest BCmetric)
optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected number of doublets (assuming ~7% doublet rate for 10X)
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.07 * ncol(pbmc))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:20, pN = 0.25, pK = optimal_pK, 
                         nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get the classification column name (it's dynamically generated)
df_col <- grep("^DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets on UMAP
DimPlot(pbmc, reduction = "umap", group.by = df_col)

# Print doublet statistics
table(pbmc@meta.data[[df_col]])

# Subset to keep only singlets
pbmc_singlets <- subset(pbmc, subset = !!sym(df_col) == "Singlet")

# Visualize singlets only
DimPlot(pbmc_singlets, reduction = "umap", group.by = "seurat_clusters")

# Print final cell count
cat("Original cells:", ncol(pbmc), "\n")
cat("Singlets retained:", ncol(pbmc_singlets), "\n")
cat("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
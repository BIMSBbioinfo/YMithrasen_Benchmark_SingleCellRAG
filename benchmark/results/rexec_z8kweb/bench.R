setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(DoubletFinder)
library(ggplot2)

# Load the PBMC dataset from working directory
# Files: barcodes.tsv, genes.tsv, matrix.mtx
sc.data <- Read10X(data.dir = ".")

# Create Seurat object with QC filters
# min.cells = 3: genes detected in at least 3 cells
# min.features = 200: cells with at least 200 detected genes
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Add mitochondrial percentage for QC
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# Filter cells based on QC metrics
# nFeature_RNA: number of genes detected per cell
# percent.mt: percentage of mitochondrial genes
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

# DoubletFinder: pK identification
sweep.res.list <- paramSweep_v3(sc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot BCmetric to identify optimal pK
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "pK Optimization", x = "pK", y = "BCmetric")

# Use optimal pK (typically the peak in BCmetric plot)
# For this example, using pK = 0.09 as commonly optimal
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Calculate expected number of doublets
# Assuming 7.5% doublet rate for ~3000 cells
annotations <- sc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
sc <- doubletFinder_v3(sc, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get the doublet classification column name (it's dynamically generated)
doublet_col <- colnames(sc@meta.data)[grep("DF.classifications", colnames(sc@meta.data))]

# Visualize doublets on UMAP
DimPlot(sc, reduction = "umap", group.by = doublet_col, pt.size = 0.5) +
  ggtitle("Doublet Detection Results") +
  theme_minimal()

# Create subset of singlets only
sc_singlets <- subset(sc, subset = sc@meta.data[[doublet_col]] == "Singlet")

# Visualize singlets only
DimPlot(sc_singlets, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5) +
  ggtitle("Singlets Only - Cluster Visualization") +
  theme_minimal()

# Summary statistics
cat("\nTotal cells:", ncol(sc), "\n")
cat("Singlets:", sum(sc@meta.data[[doublet_col]] == "Singlet"), "\n")
cat("Doublets:", sum(sc@meta.data[[doublet_col]] == "Doublet"), "\n")
cat("Doublet rate:", round(sum(sc@meta.data[[doublet_col]] == "Doublet") / ncol(sc) * 100, 2), "%\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
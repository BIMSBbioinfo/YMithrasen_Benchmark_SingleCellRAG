setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(DoubletFinder)

# Load 10x data
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

# Add mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filter cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Standard preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)

# DoubletFinder parameter sweep
sweep.res <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Get optimal pK
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Estimate doublet number (assuming 7.5% doublet rate for ~10k cells recovered)
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(pbmc@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, PCs = 1:30, pN = 0.25, pK = pK_optimal, 
                         nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

# Get classification column name
df_col <- colnames(pbmc@meta.data)[grep("DF.classifications", colnames(pbmc@meta.data))]

# Visualize doublets
DimPlot(pbmc, reduction = "umap", group.by = df_col)

# Remove doublets
pbmc_singlets <- subset(pbmc, subset = pbmc@meta.data[[df_col]] == "Singlet")

# Verify singlet count
print(paste("Original cells:", ncol(pbmc)))
print(paste("Singlets retained:", ncol(pbmc_singlets)))
print(paste("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets)))

# Visualize singlets only
DimPlot(pbmc_singlets, reduction = "umap", group.by = "seurat_clusters")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Load PBMC data (assuming 10X Genomics format)
pbmc.data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC", min.cells = 3, min.features = 200)

# 2. Add metadata: percent mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 3. QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 4. Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 5. Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 6. Scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 7. PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 8. Determine dimensionality
ElbowPlot(pbmc)

# 9. UMAP
pbmc <- RunUMAP(pbmc, dims = 1:20)

# 10. Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 11. Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 12. Identify cluster markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 13. Save markers
write.csv(pbmc.markers, file = "PBMC_cluster_markers.csv")

# 14. Annotate clusters using known markers
cluster.ids <- c("CD4 T", "CD8 T", "NK", "B", "Monocytes", "Dendritic", "Megakaryocytes")
names(cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cluster.ids)

# 15. Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 16. Save Seurat object
saveRDS(pbmc, file = "PBMC_seurat_object.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
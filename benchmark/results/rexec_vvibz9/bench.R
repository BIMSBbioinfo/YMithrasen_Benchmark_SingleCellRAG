setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)

# Load the data
pbmc.data <- Read10X(data.dir = ".")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)

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

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Find marker genes for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")

# Get top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top_markers)

# Annotate cell types using SingleR
ref <- celldex::MonacoImmuneData()
pbmc_counts <- GetAssayData(pbmc, slot = "data")
singler_results <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)

# Add SingleR annotations to Seurat object
pbmc$singler_labels <- singler_results$labels

# Visualize clusters with original identities
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + 
  ggtitle("PBMC Clusters")

# Visualize with SingleR annotations
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "singler_labels", label = TRUE, label.size = 3) + 
  ggtitle("SingleR Cell Type Annotations")

print(p1)
print(p2)

# Save results
write.csv(pbmc.markers, "pbmc_cluster_markers.csv", row.names = FALSE)
write.csv(data.frame(cell = colnames(pbmc), cluster = Idents(pbmc), cell_type = pbmc$singler_labels), 
          "pbmc_cell_annotations.csv", row.names = FALSE)

saveRDS(pbmc, "pbmc_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
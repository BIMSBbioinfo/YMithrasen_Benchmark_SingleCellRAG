setwd("D:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(Matrix)
library(ggplot2)

# Step 1: Load the data
# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes and genes
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes$V2  # Use HGNC symbols
colnames(matrix) <- barcodes$V1

# Step 2: Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

# Step 3: Preprocessing
# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Step 4: Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), 
               npcs = 50, verbose = FALSE)

# Visualize PCA results (optional)
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Step 5: Build KNN graph using 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20)

# Step 6: Cluster cells with resolution 0.5
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1)

# Step 7: Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Step 8: Visualize clusters on UMAP
umap_plot <- DimPlot(pbmc, reduction = "umap", 
                     label = TRUE, 
                     pt.size = 0.5,
                     label.size = 6) +
  ggtitle("PBMC Clusters (2,700 cells)") +
  theme_minimal()

print(umap_plot)

# Save the plot
ggsave("pbmc_umap_clusters.png", plot = umap_plot, 
       width = 10, height = 8, dpi = 300)

# Print cluster summary
cat("\nNumber of cells per cluster:\n")
print(table(Idents(pbmc)))

# Optional: Save the Seurat object for future analysis
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

# Additional visualization: UMAP colored by number of genes detected
FeaturePlot(pbmc, features = "nFeature_RNA", reduction = "umap")

# Show cluster markers (top 2 genes per cluster)
cat("\nFinding cluster markers...\n")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25,
                               max.cells.per.ident = 100)

# Display top markers
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
print(top_markers)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
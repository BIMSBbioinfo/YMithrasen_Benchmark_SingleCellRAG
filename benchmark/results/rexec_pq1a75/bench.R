setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(Matrix)

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, project = "PBMC", genes = genes[, 2], cells = barcodes[, 1])

# Filter cells and genes
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Scale data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine significant PCs
elbow <- InflectionPoint(pbmc@reductions$pca@stddev)
pbmc <- RunUMAP(pbmc, dims = 1:elbow)

# Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:elbow)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Perform differential expression analysis
cluster0 <- subset(pbmc, idents = "0")
cluster3 <- subset(pbmc, idents = "3")
diff <- FindMarkers(cluster0, cluster3, only.pos = FALSE)

# Create volcano plot
volcano_data <- data.frame(log2FC = diff$avg_logFC, p_val_adj = diff$p_val_adj)
volcano_data$significant <- ifelse(abs(volcano_data$log2FC) > 1 & volcano_data$p_val_adj < 0.05, "yes", "no")

ggplot(volcano_data, aes(x = log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
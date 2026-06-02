setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load necessary libraries
library(Seurat)
library(Matrix)
library(ggplot2)

# Load data
barcodes <- read.table("barcodes.tsv", colClasses = "character")
genes <- read.table("genes.tsv", colClasses = c("character", "character"))
matrix <- readMM("matrix.mtx")

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = matrix, 
                                    project = "PBMCs", 
                                    genes = genes[, 2], 
                                    cells = barcodes[, 1])

# Normalize data
seurat_object <- NormalizeData(seurat_object, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# Find variable features
seurat_object <- FindVariableFeatures(seurat_object, 
                                      selection.method = "vst")

# Scale data
seurat_object <- ScaleData(seurat_object, 
                            features = VariableFeatures(seurat_object))

# Run PCA
seurat_object <- RunPCA(seurat_object, 
                         features = VariableFeatures(seurat_object))

# Determine significant PCs
seurat_object <- JackStraw(seurat_object, 
                          num.replicate = 100)

# Select significant PCs
seurat_object <- RunUMAP(seurat_object, 
                         dims = 1:10)

# Cluster cells
seurat_object <- FindNeighbors(seurat_object, 
                                dims = 1:10)
seurat_object <- FindClusters(seurat_object, 
                              resolution = 0.5)

# Perform differential expression analysis
cluster0 <- subset(seurat_object, idents = "0")
cluster3 <- subset(seurat_object, idents = "3")
diff_exp <- FindMarkers(cluster0, cluster3)

# Create volcano plot
volcano_data <- data.frame(log2FC = diff_exp$avg_logFC, 
                           p_val_adj = diff_exp$p_val_adj)
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
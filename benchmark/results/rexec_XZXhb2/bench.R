setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("Matrix", quietly=TRUE)) tryCatch(BiocManager::install("Matrix", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Matrix", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Load the 10X data
matrix <- readMM("matrix.mtx")
barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
genes <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)

# Set row and column names
rownames(matrix) <- genes[, 2]
colnames(matrix) <- barcodes[, 1]

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, project = "PBMC", min.cells = 3, min.features = 200)

# Normalize data
pbmc <- NormalizeData(pbmc)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Verify clusters 0 and 3 exist
if(!all(c("0", "3") %in% levels(Idents(pbmc)))) {
  stop("Clusters 0 and 3 not found. Adjust resolution parameter.")
}

# Differential expression analysis between cluster 0 and cluster 3
de_results <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 3)

# Add gene names as a column
de_results$gene <- rownames(de_results)

# Remove rows with NA or Inf values in log2FC or p_val_adj
de_results <- de_results[!is.na(de_results$avg_log2FC) & 
                         !is.infinite(de_results$avg_log2FC) & 
                         !is.na(de_results$p_val_adj) & 
                         !is.infinite(de_results$p_val_adj), ]

# Add significance column
de_results$significant <- ifelse(abs(de_results$avg_log2FC) > 1 & 
                                 de_results$p_val_adj < 0.05, 
                                 "Significant", "Not Significant")

# Create volcano plot
volcano_plot <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value",
       title = "Volcano Plot: Cluster 0 vs Cluster 3",
       color = "Significance") +
  theme_minimal() +
  theme(legend.position = "right")

print(volcano_plot)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
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
library(dplyr)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data
# ============================================================================

# Read the Matrix Market format file
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# STEP 2: Quality Control
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 3: Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# STEP 4: Feature Selection
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================================================
# STEP 5: Scaling
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 6: PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# STEP 7: Clustering
# ============================================================================

# Find neighbors
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# STEP 8: UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

# ============================================================================
# STEP 9: Differential Expression Analysis (Cluster 0 vs Cluster 3)
# ============================================================================

# Check if clusters 0 and 3 exist
if (!("0" %in% levels(Idents(pbmc))) || !("3" %in% levels(Idents(pbmc)))) {
  stop("Clusters 0 or 3 not found. Available clusters: ", 
       paste(levels(Idents(pbmc)), collapse = ", "))
}

# Perform differential expression
de_results <- FindMarkers(pbmc, 
                          ident.1 = "0", 
                          ident.2 = "3", 
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          test.use = "wilcox")

# Add gene names as a column
de_results$gene <- rownames(de_results)

# Calculate -log10(adjusted p-value)
de_results$neg_log10_padj <- -log10(de_results$p_val_adj)

# Replace infinite values (from p_val_adj = 0) with a large number
de_results$neg_log10_padj[is.infinite(de_results$neg_log10_padj)] <- 
  max(de_results$neg_log10_padj[!is.infinite(de_results$neg_log10_padj)]) + 10

# Classify genes as significant or not
de_results$significant <- ifelse(abs(de_results$avg_log2FC) > 1 & 
                                  de_results$p_val_adj < 0.05, 
                                  "Significant", 
                                  "Not Significant")

# Count significant genes
n_sig <- sum(de_results$significant == "Significant")
cat("\nNumber of significantly differentially expressed genes:", n_sig, "\n")
cat("(|log2FC| > 1 and adjusted p-value < 0.05)\n")

# Display top differentially expressed genes
cat("\nTop 10 upregulated genes in cluster 0:\n")
print(head(de_results[order(-de_results$avg_log2FC), ], 10))

cat("\nTop 10 downregulated genes in cluster 0:\n")
print(head(de_results[order(de_results$avg_log2FC), ], 10))

# ============================================================================
# STEP 10: Volcano Plot
# ============================================================================

# Create volcano plot
volcano_plot <- ggplot(de_results, aes(x = avg_log2FC, 
                                       y = neg_log10_padj, 
                                       color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", 
                                 "Significant" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Differential Expression: Cluster 0 vs Cluster 3",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

# Display the plot
print(volcano_plot)

# Save the plot
ggsave("volcano_plot_cluster0_vs_cluster3.pdf", 
       plot = volcano_plot, 
       width = 10, 
       height = 8)

ggsave("volcano_plot_cluster0_vs_cluster3.png", 
       plot = volcano_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

# ============================================================================
# STEP 11: Save Results
# ============================================================================

# Save differential expression results
write.csv(de_results, 
          "DE_cluster0_vs_cluster3.csv", 
          row.names = FALSE)

# Save Seurat object
saveRDS(pbmc, file = "pbmc_seurat_object.rds")

cat("\nAnalysis complete!\n")
cat("Results saved:\n")
cat("  - Volcano plot: volcano_plot_cluster0_vs_cluster3.pdf/png\n")
cat("  - DE results: DE_cluster0_vs_cluster3.csv\n")
cat("  - Seurat object: pbmc_seurat_object.rds\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
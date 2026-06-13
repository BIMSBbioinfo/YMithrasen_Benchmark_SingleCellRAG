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
library(Matrix)

# ============================================================================
# 1. LOAD DATA
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
# 2. QUALITY CONTROL
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("QC_violin_plots.pdf", width = 12, height = 4)

# Feature scatter plots
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("QC_scatter_plots.pdf", width = 10, height = 4)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# 3. NORMALIZATION
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# 4. FEATURE SELECTION
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
ggsave("variable_features.pdf", width = 10, height = 6)

# ============================================================================
# 5. SCALING
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# 6. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave("PCA_loadings.pdf", width = 10, height = 8)

DimPlot(pbmc, reduction = "pca")
ggsave("PCA_plot.pdf", width = 8, height = 6)

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)
ggsave("elbow_plot.pdf", width = 8, height = 6)

# ============================================================================
# 7. CLUSTERING
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Check cluster assignments
cat("Number of clusters:", length(unique(Idents(pbmc))), "\n")
table(Idents(pbmc))

# ============================================================================
# 8. NON-LINEAR DIMENSIONALITY REDUCTION (UMAP)
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("UMAP_clusters.pdf", width = 8, height = 6)

# ============================================================================
# 9. DIFFERENTIAL EXPRESSION: CLUSTER 0 vs CLUSTER 3
# ============================================================================

# Find markers between cluster 0 and cluster 3
cluster0_vs_3 <- FindMarkers(pbmc, 
                             ident.1 = 0, 
                             ident.2 = 3, 
                             min.pct = 0.25,
                             logfc.threshold = 0.25)

# Add gene names as a column
cluster0_vs_3$gene <- rownames(cluster0_vs_3)

# View top differentially expressed genes
cat("\nTop 10 upregulated genes in cluster 0 vs 3:\n")
print(head(cluster0_vs_3[order(cluster0_vs_3$avg_log2FC, decreasing = TRUE), ], 10))

cat("\nTop 10 downregulated genes in cluster 0 vs 3:\n")
print(head(cluster0_vs_3[order(cluster0_vs_3$avg_log2FC, decreasing = FALSE), ], 10))

# Save results
write.csv(cluster0_vs_3, "cluster0_vs_3_DEGs.csv", row.names = FALSE)

# ============================================================================
# 10. VOLCANO PLOT
# ============================================================================

# Prepare data for volcano plot
cluster0_vs_3$diffexpressed <- "NO"
cluster0_vs_3$diffexpressed[cluster0_vs_3$avg_log2FC > 1 & cluster0_vs_3$p_val_adj < 0.05] <- "UP"
cluster0_vs_3$diffexpressed[cluster0_vs_3$avg_log2FC < -1 & cluster0_vs_3$p_val_adj < 0.05] <- "DOWN"

# Add labels for top genes
cluster0_vs_3$delabel <- NA
cluster0_vs_3$delabel[cluster0_vs_3$diffexpressed != "NO"] <- 
  cluster0_vs_3$gene[cluster0_vs_3$diffexpressed != "NO"]

# Select top 10 genes by p-value for labeling
top_genes <- cluster0_vs_3 %>%
  filter(diffexpressed != "NO") %>%
  arrange(p_val_adj) %>%
  head(10)

cluster0_vs_3$delabel <- ifelse(cluster0_vs_3$gene %in% top_genes$gene, 
                                 cluster0_vs_3$gene, NA)

# Create volcano plot
volcano_plot <- ggplot(data = cluster0_vs_3, 
                       aes(x = avg_log2FC, 
                           y = -log10(p_val_adj), 
                           col = diffexpressed,
                           label = delabel)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                     labels = c("DOWN" = "Downregulated", 
                               "NO" = "Not significant", 
                               "UP" = "Upregulated")) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = "Differential Expression: Cluster 0 vs Cluster 3",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Expression") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right") +
  ggrepel::geom_text_repel(max.overlaps = 15, 
                          size = 3,
                          box.padding = 0.5,
                          point.padding = 0.3)

print(volcano_plot)
ggsave("volcano_plot_cluster0_vs_3.pdf", width = 10, height = 8)

# ============================================================================
# 11. SUMMARY STATISTICS
# ============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
cat("Total genes tested:", nrow(cluster0_vs_3), "\n")
cat("Significantly upregulated (|log2FC| > 1, p_adj < 0.05):", 
    sum(cluster0_vs_3$diffexpressed == "UP"), "\n")
cat("Significantly downregulated (|log2FC| > 1, p_adj < 0.05):", 
    sum(cluster0_vs_3$diffexpressed == "DOWN"), "\n")
cat("Not significant:", sum(cluster0_vs_3$diffexpressed == "NO"), "\n")

# ============================================================================
# 12. SAVE SEURAT OBJECT
# ============================================================================

saveRDS(pbmc, file = "pbmc_seurat_object.rds")
cat("\nSeurat object saved as 'pbmc_seurat_object.rds'\n")

cat("\nAnalysis complete!\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
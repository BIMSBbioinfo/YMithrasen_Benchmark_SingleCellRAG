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
library(EnhancedVolcano)
library(pheatmap)

# Assuming you have a Seurat object named 'pbmc'
# If not, load your data first:
# pbmc <- readRDS("path/to/your/seurat_object.rds")

# ============================================
# 1. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================

# Define clusters to compare (e.g., cluster 0 vs cluster 1)
cluster1 <- 0
cluster2 <- 1

# Perform differential expression using Wilcoxon rank sum test
de_markers <- FindMarkers(
  pbmc,
  ident.1 = cluster1,
  ident.2 = cluster2,
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = FALSE
)

# Add gene names as a column
de_markers$gene <- rownames(de_markers)

# Filter significant genes (adjusted p-value < 0.05)
sig_genes <- de_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

# View top upregulated and downregulated genes
head(sig_genes, 20)

# ============================================
# 2. VOLCANO PLOT
# ============================================

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  de_markers,
  lab = rownames(de_markers),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = paste0('Cluster ', cluster1, ' vs Cluster ', cluster2),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 2.0,
  labSize = 4.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.5,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0
)

print(volcano_plot)
ggsave("volcano_plot.pdf", volcano_plot, width = 10, height = 8)

# ============================================
# 3. HEATMAP OF TOP GENES
# ============================================

# Select top 20 upregulated and downregulated genes
top_genes <- sig_genes %>%
  group_by(avg_log2FC > 0) %>%
  slice_max(order_by = abs(avg_log2FC), n = 20) %>%
  pull(gene)

# Extract expression data for heatmap
heatmap_data <- AverageExpression(
  pbmc,
  features = top_genes,
  group.by = "seurat_clusters"
)$RNA

# Create heatmap
pdf("de_heatmap.pdf", width = 8, height = 10)
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top DE Genes Heatmap",
  fontsize = 10
)
dev.off()

# ============================================
# 4. VIOLIN PLOTS FOR TOP GENES
# ============================================

# Select top 6 genes for visualization
top6_genes <- head(sig_genes$gene, 6)

# Create violin plots
violin_plots <- VlnPlot(
  pbmc,
  features = top6_genes,
  idents = c(cluster1, cluster2),
  ncol = 3,
  pt.size = 0.1
)

ggsave("violin_plots.pdf", violin_plots, width = 12, height = 8)

# ============================================
# 5. FEATURE PLOTS (UMAP)
# ============================================

# Visualize top genes on UMAP
feature_plots <- FeaturePlot(
  pbmc,
  features = top6_genes,
  ncol = 3,
  reduction = "umap"
)

ggsave("feature_plots.pdf", feature_plots, width = 12, height = 8)

# ============================================
# 6. DOT PLOT
# ============================================

# Create dot plot for top genes
dot_plot <- DotPlot(
  pbmc,
  features = top6_genes,
  idents = c(cluster1, cluster2)
) + 
  RotatedAxis() +
  ggtitle("Expression of Top DE Genes")

ggsave("dot_plot.pdf", dot_plot, width = 10, height = 6)

# ============================================
# 7. EXPORT RESULTS
# ============================================

# Save all DE results
write.csv(de_markers, "de_results_all.csv", row.names = FALSE)

# Save significant genes only
write.csv(sig_genes, "de_results_significant.csv", row.names = FALSE)

# Summary statistics
cat("\n=== Differential Expression Summary ===\n")
cat("Total genes tested:", nrow(de_markers), "\n")
cat("Significant genes (p_adj < 0.05):", nrow(sig_genes), "\n")
cat("Upregulated in cluster", cluster1, ":", 
    sum(sig_genes$avg_log2FC > 0), "\n")
cat("Downregulated in cluster", cluster1, ":", 
    sum(sig_genes$avg_log2FC < 0), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
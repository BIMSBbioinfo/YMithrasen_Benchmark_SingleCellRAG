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
# If loading from file:
# pbmc <- readRDS("pbmc_seurat_object.rds")

# ============================================
# 1. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================

# Define clusters to compare (e.g., cluster 0 vs cluster 1)
cluster1 <- 0
cluster2 <- 1

# Perform differential expression using Wilcoxon Rank Sum test
de_markers <- FindMarkers(
  object = pbmc,
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
sig_markers <- de_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

# View top upregulated and downregulated genes
top_up <- head(sig_markers %>% filter(avg_log2FC > 0), 10)
top_down <- head(sig_markers %>% filter(avg_log2FC < 0), 10)

print("Top 10 Upregulated Genes:")
print(top_up)
print("Top 10 Downregulated Genes:")
print(top_down)

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
ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)

# ============================================
# 3. HEATMAP OF TOP MARKERS
# ============================================

# Select top 20 differentially expressed genes
top_genes <- sig_markers %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20) %>%
  pull(gene)

# Extract expression data for these genes
Idents(pbmc) <- pbmc$seurat_clusters
cells_cluster1 <- WhichCells(pbmc, idents = cluster1)
cells_cluster2 <- WhichCells(pbmc, idents = cluster2)

# Get average expression per cluster
avg_exp <- AverageExpression(
  pbmc,
  features = top_genes,
  group.by = "seurat_clusters"
)$RNA

# Create heatmap
pheatmap(
  avg_exp[, c(as.character(cluster1), as.character(cluster2))],
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 20 DE Genes",
  fontsize = 10,
  filename = "heatmap_top_genes.png",
  width = 6,
  height = 10
)

# ============================================
# 4. VIOLIN PLOTS FOR TOP GENES
# ============================================

# Create violin plots for top 6 genes
top_6_genes <- head(top_genes, 6)

violin_plots <- VlnPlot(
  pbmc,
  features = top_6_genes,
  idents = c(cluster1, cluster2),
  ncol = 3,
  pt.size = 0.1
)

print(violin_plots)
ggsave("violin_plots_top_genes.png", violin_plots, width = 12, height = 8, dpi = 300)

# ============================================
# 5. FEATURE PLOTS (UMAP)
# ============================================

feature_plot <- FeaturePlot(
  pbmc,
  features = top_6_genes,
  ncol = 3,
  reduction = "umap"
)

print(feature_plot)
ggsave("feature_plots_umap.png", feature_plot, width = 12, height = 8, dpi = 300)

# ============================================
# 6. DOT PLOT
# ============================================

dot_plot <- DotPlot(
  pbmc,
  features = top_genes,
  idents = c(cluster1, cluster2),
  cols = c("lightgrey", "red"),
  dot.scale = 8
) + RotatedAxis()

print(dot_plot)
ggsave("dot_plot_top_genes.png", dot_plot, width = 10, height = 8, dpi = 300)

# ============================================
# 7. EXPORT RESULTS
# ============================================

# Save all DE results
write.csv(de_markers, "DE_results_all.csv", row.names = FALSE)
write.csv(sig_markers, "DE_results_significant.csv", row.names = FALSE)

# Create summary statistics
summary_stats <- data.frame(
  Total_Genes_Tested = nrow(de_markers),
  Significant_Genes = nrow(sig_markers),
  Upregulated = sum(sig_markers$avg_log2FC > 0),
  Downregulated = sum(sig_markers$avg_log2FC < 0),
  Mean_LogFC = mean(sig_markers$avg_log2FC)
)

write.csv(summary_stats, "DE_summary_statistics.csv", row.names = FALSE)

print("Analysis complete! All plots and results saved.")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
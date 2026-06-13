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
top_up <- head(sig_genes %>% filter(avg_log2FC > 0), 10)
top_down <- head(sig_genes %>% filter(avg_log2FC < 0), 10)

print("Top 10 upregulated genes:")
print(top_up)
print("Top 10 downregulated genes:")
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
# 3. HEATMAP OF TOP GENES
# ============================================

# Select top 20 differentially expressed genes
top_genes <- head(sig_genes, 20)$gene

# Extract expression data for these genes
expr_data <- GetAssayData(pbmc, slot = "scale.data")[top_genes, ]

# Get cell identities
cell_clusters <- pbmc@active.ident

# Filter for cells in the two clusters of interest
cells_of_interest <- names(cell_clusters)[cell_clusters %in% c(cluster1, cluster2)]
expr_subset <- expr_data[, cells_of_interest]

# Create annotation for heatmap
annotation_col <- data.frame(
  Cluster = cell_clusters[cells_of_interest]
)
rownames(annotation_col) <- cells_of_interest

# Plot heatmap
pheatmap(
  expr_subset,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = paste0("Top 20 DE genes: Cluster ", cluster1, " vs ", cluster2),
  fontsize_row = 8,
  filename = "de_heatmap.png",
  width = 10,
  height = 8
)

# ============================================
# 4. VIOLIN PLOTS FOR TOP GENES
# ============================================

# Create violin plots for top 6 genes
top_6_genes <- head(sig_genes, 6)$gene

violin_plot <- VlnPlot(
  pbmc,
  features = top_6_genes,
  idents = c(cluster1, cluster2),
  ncol = 3,
  pt.size = 0.1
)

print(violin_plot)
ggsave("violin_plots.png", violin_plot, width = 12, height = 8, dpi = 300)

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
ggsave("feature_plots.png", feature_plot, width = 12, height = 8, dpi = 300)

# ============================================
# 6. DOT PLOT
# ============================================

dot_plot <- DotPlot(
  pbmc,
  features = top_genes[1:15],
  idents = c(cluster1, cluster2)
) + 
  RotatedAxis() +
  ggtitle(paste0("Top 15 DE genes: Cluster ", cluster1, " vs ", cluster2))

print(dot_plot)
ggsave("dot_plot.png", dot_plot, width = 10, height = 6, dpi = 300)

# ============================================
# 7. EXPORT RESULTS
# ============================================

# Save all DE results to CSV
write.csv(de_markers, "de_results_all.csv", row.names = TRUE)
write.csv(sig_genes, "de_results_significant.csv", row.names = TRUE)

# Create summary statistics
summary_stats <- data.frame(
  Total_genes_tested = nrow(de_markers),
  Significant_genes = nrow(sig_genes),
  Upregulated = sum(sig_genes$avg_log2FC > 0),
  Downregulated = sum(sig_genes$avg_log2FC < 0)
)

write.csv(summary_stats, "de_summary.csv", row.names = FALSE)

print("Analysis complete! Results saved.")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
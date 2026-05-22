setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))

# Set default assay to RNA
DefaultAssay(pbmc_clustered) <- "RNA"

# Verify clusters exist
if(!all(c("0", "3") %in% levels(Idents(pbmc_clustered)))) {
  stop("Clusters 0 and/or 3 not found in the Seurat object")
}

# Perform differential expression analysis between cluster 0 and cluster 3
cluster0_vs_3_markers <- FindMarkers(
  pbmc_clustered,
  ident.1 = 0,
  ident.2 = 3,
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0
)

# Add gene names as a column
cluster0_vs_3_markers$gene <- rownames(cluster0_vs_3_markers)

# Filter out genes with p_val_adj = 0 to avoid Inf values
cluster0_vs_3_markers <- cluster0_vs_3_markers %>%
  filter(p_val_adj > 0)

# Add significance category for coloring
cluster0_vs_3_markers <- cluster0_vs_3_markers %>%
  mutate(
    significance = case_when(
      abs(avg_log2FC) > 1 & p_val_adj < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Create volcano plot
volcano_plot <- ggplot(cluster0_vs_3_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Cluster 0 vs Cluster 3",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(volcano_plot)

# Display summary of significant genes
sig_genes <- cluster0_vs_3_markers %>%
  filter(significance == "Significant") %>%
  arrange(p_val_adj)

cat("\nNumber of significant genes (|log2FC| > 1 & p_val_adj < 0.05):", nrow(sig_genes), "\n")
cat("\nTop 10 significant genes:\n")
print(head(sig_genes[, c("gene", "avg_log2FC", "p_val_adj")], 10))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
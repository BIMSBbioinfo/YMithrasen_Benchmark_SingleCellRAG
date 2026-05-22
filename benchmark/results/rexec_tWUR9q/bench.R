setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Perform differential expression between cluster 0 and cluster 3
deg_results <- FindMarkers(pbmc_clustered, 
                          ident.1 = 0, 
                          ident.2 = 3, 
                          test.use = "wilcox")

# Add gene names as a column
deg_results$gene <- rownames(deg_results)

# Calculate -log10(adjusted p-value)
deg_results$neg_log10_padj <- -log10(deg_results$p_val_adj)

# Add color column based on thresholds
deg_results$color <- ifelse(abs(deg_results$avg_log2FC) > 1 & deg_results$p_val_adj < 0.05, 
                            "red", 
                            "grey")

# Create volcano plot
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = color), alpha = 0.6) +
  scale_color_identity() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "log2 Fold Change", 
       y = "-log10(Adjusted P-value)",
       title = "Volcano Plot: Cluster 0 vs Cluster 3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(volcano_plot)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
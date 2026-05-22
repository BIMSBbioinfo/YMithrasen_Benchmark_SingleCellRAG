setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Perform differential expression analysis between cluster 0 and cluster 3
deg_results <- FindMarkers(pbmc_clustered, 
                          ident.1 = 0, 
                          ident.2 = 3, 
                          test.use = "wilcox",
                          logfc.threshold = 0)

# Add gene names from rownames
deg_results$gene <- rownames(deg_results)

# Filter out NA values in p_val_adj
deg_results <- deg_results[!is.na(deg_results$p_val_adj), ]

# Replace p_val_adj = 0 with smallest non-zero value to avoid Inf
min_pval <- min(deg_results$p_val_adj[deg_results$p_val_adj > 0])
deg_results$p_val_adj[deg_results$p_val_adj == 0] <- min_pval

# Calculate -log10(p_val_adj)
deg_results$neg_log10_padj <- -log10(deg_results$p_val_adj)

# Create significance category
deg_results$significance <- ifelse(abs(deg_results$avg_log2FC) > 1 & deg_results$p_val_adj < 0.05, 
                                   "Significant", 
                                   "Not Significant")

# Create volcano plot
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_padj, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "log2 Fold Change", 
       y = "-log10(adjusted p-value)",
       title = "Volcano Plot: Cluster 0 vs Cluster 3") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(volcano_plot)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
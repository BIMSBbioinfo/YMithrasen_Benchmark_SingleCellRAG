setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))

# Visualize marker genes to identify cluster identities
FeaturePlot(pbmc, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7"), ncol = 4)

# Violin plots for marker gene expression across clusters
VlnPlot(pbmc, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7"), ncol = 4)

# Based on marker expression, create cell type labels
# Note: Adjust cluster numbers based on your actual marker expression patterns
new.cluster.ids <- c("CD4 T cells", "CD14 Monocytes", "CD4 T cells", "B cells", 
                     "CD8 T cells", "NK cells", "CD14 Monocytes", "NK cells")
names(new.cluster.ids) <- levels(pbmc)

# Rename cluster identities to cell types
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Generate UMAP plot with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE) + 
  ggtitle("PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
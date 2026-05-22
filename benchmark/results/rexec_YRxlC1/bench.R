setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Visualize canonical markers to identify cell types
markers_to_check <- c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7")

# Feature plot to see marker expression across clusters
FeaturePlot(pbmc, features = markers_to_check, ncol = 4)

# Violin plot for clearer cluster-specific expression
VlnPlot(pbmc, features = markers_to_check, ncol = 4)

# Dot plot for comprehensive view
DotPlot(pbmc, features = markers_to_check) + RotatedAxis()

# After inspecting the plots above, create mapping
# NOTE: Adjust cluster numbers based on YOUR marker expression patterns
# This is an example - verify with your plots above
new.cluster.ids <- c(
  "Naive CD4 T",     # Cluster with high IL7R, CCR7
  "CD14+ Mono",      # Cluster with high CD14, LYZ
  "Memory CD4 T",    # Cluster with IL7R but lower CCR7
  "B",               # Cluster with high MS4A1
  "CD8 T",           # Cluster with high CD8A
  "FCGR3A+ Mono",    # Another monocyte cluster
  "NK",              # Cluster with high GNLY, NKG7
  "DC",              # Dendritic cells if present
  "Platelet"         # Platelets if present
)

# Get current cluster levels
current_levels <- levels(pbmc)
print(paste("Current cluster levels:", paste(current_levels, collapse = ", ")))

# Assign names to match current levels
names(new.cluster.ids) <- levels(pbmc)

# Rename cluster identities
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Store cell type annotations in metadata
pbmc$cell_type <- Idents(pbmc)

# Generate UMAP with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5, repel = TRUE) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  ggtitle("PBMC Cell Types") +
  theme(
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# Alternative: UMAP without cluster numbers
DimPlot(pbmc, reduction = "umap", label = FALSE) + 
  ggtitle("PBMC Cell Types - Labeled")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
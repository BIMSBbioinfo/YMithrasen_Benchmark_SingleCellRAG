setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
library(Seurat)
library(dplyr)
library(ggplot2)

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, 
                               test.use = "wilcox",
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# Filter by adjusted p-value
pbmc.markers <- pbmc.markers %>% 
  filter(p_val_adj < 0.05)

# View top markers per cluster
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Assign cell type identities based on canonical markers
# CD3D/CD3E: T cells
# CD14/LYZ: CD14+ Monocytes
# MS4A1/CD79A: B cells
# NKG7/GNLY: NK cells
# FCGR3A: CD16+ Monocytes
# IL7R: Naive CD4+ T
# CCR7: Memory CD4+ T
# CD8A: CD8+ T cells

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                     "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Visualize annotated clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))

# Visualize key marker genes
FeaturePlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7", 
                                "FCGR3A", "IL7R", "CD8A"))

# Violin plot of marker genes
VlnPlot(pbmc, features = c("CD3D", "CD14", "MS4A1", "NKG7", 
                           "FCGR3A", "IL7R", "CD8A"))

# Heatmap of top markers
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))

# Load the Seurat object (assuming it's saved as an RDS file)
# If your object is already in memory with a different name, skip this line
pbmc <- readRDS("pbmc_seurat_object.rds")

# Visualize marker genes to confirm cluster identities
FeaturePlot(pbmc, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7"), ncol = 4)

# Check current cluster levels
levels(pbmc)

# Create named vector mapping cluster numbers to cell types
# Adjust the mapping based on your actual cluster-to-celltype correspondence
# This is an example - you need to verify which cluster number corresponds to which cell type
new.cluster.ids <- c("CD4 T cells", "CD14 Monocytes", "CD4 T cells", "B cells", 
                     "CD8 T cells", "NK cells", "CD14 Monocytes", "CD8 T cells", 
                     "B cells")
names(new.cluster.ids) <- levels(pbmc)

# Rename cluster identities
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Store cell type identities in metadata
pbmc$celltype <- Idents(pbmc)

# Reorder levels for better visualization
pbmc <- SetIdent(pbmc, value = "celltype")
levels(pbmc) <- c("CD4 T cells", "CD8 T cells", "B cells", "CD14 Monocytes", "NK cells")

# Generate UMAP plot with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(Seurat)))

# Load the Seurat object (assuming it's saved as an RDS file)
# If your object is already in memory with a different name, skip this line
pbmc <- readRDS("pbmc_seurat_object.rds")

# Verify current cluster identities
table(Idents(pbmc))

# Define new cluster identities based on marker gene expression
# Adjust cluster numbers based on your actual marker expression patterns
new.cluster.ids <- c("CD4 T cells", "CD14 Monocytes", "CD4 T cells", "B cells", 
                     "CD8 T cells", "CD14 Monocytes", "NK cells", "CD4 T cells", 
                     "B cells")

# Get current cluster levels
names(new.cluster.ids) <- levels(Idents(pbmc))

# Rename cluster identities
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Generate UMAP plot with cell type labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend()

# Alternative plot with legend
DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5, pt.size = 0.5)

# Save the annotated Seurat object
saveRDS(pbmc, "pbmc_annotated.rds")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
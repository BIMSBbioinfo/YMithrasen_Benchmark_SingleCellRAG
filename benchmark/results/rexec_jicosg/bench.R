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
library(Matrix)

# Set working directory (adjust as needed)
# setwd("path/to/your/data")

# ============================================================================
# STEP 1: Load the data
# ============================================================================

# Read the Matrix Market format sparse matrix
matrix <- readMM("matrix.mtx")

# Read barcodes (cell identifiers)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene identifiers and symbols)
genes <- read.table("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("gene_id", "gene_symbol")

# Assign row and column names to the matrix
rownames(matrix) <- genes$gene_symbol
colnames(matrix) <- barcodes$barcode

# Create Seurat object
pbmc <- CreateSeuratObject(counts = matrix, 
                           project = "PBMC2700", 
                           min.cells = 3, 
                           min.features = 200)

cat("Loaded data:", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

# ============================================================================
# STEP 2: Quality Control
# ============================================================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("After QC filtering:", ncol(pbmc), "cells remaining\n")

# ============================================================================
# STEP 3: Normalization
# ============================================================================

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================================
# STEP 4: Feature Selection
# ============================================================================

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# ============================================================================
# STEP 5: Scaling
# ============================================================================

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ============================================================================
# STEP 6: PCA
# ============================================================================

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualize PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================================================
# STEP 7: Clustering
# ============================================================================

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

cat("Number of clusters identified:", length(unique(Idents(pbmc))), "\n")

# ============================================================================
# STEP 8: UMAP
# ============================================================================

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("PBMC Clusters")

# ============================================================================
# STEP 9: Find Cluster Markers
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print(top_markers)

# ============================================================================
# STEP 10: Cell Type Annotation Based on Canonical Markers
# ============================================================================

# Define canonical marker genes for PBMC cell types
canonical_markers <- list(
  "CD4 T cells" = c("IL7R", "CCR7"),
  "CD14 Monocytes" = c("CD14", "LYZ"),
  "B cells" = c("MS4A1"),
  "CD8 T cells" = c("CD8A"),
  "NK cells" = c("GNLY", "NKG7")
)

# Visualize canonical markers
FeaturePlot(pbmc, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", 
                                "CD8A", "GNLY", "NKG7"), 
            ncol = 4)

# Dot plot of canonical markers by cluster
DotPlot(pbmc, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", 
                           "CD8A", "GNLY", "NKG7")) + 
  RotatedAxis() +
  ggtitle("Canonical Marker Expression by Cluster")

# Calculate average expression of markers per cluster
avg_exp <- AverageExpression(pbmc, features = unlist(canonical_markers), return.seurat = FALSE)
print(avg_exp$RNA)

# ============================================================================
# STEP 11: Assign Cell Type Labels
# ============================================================================

# Manual annotation based on marker expression
# Adjust cluster numbers based on your actual results
new.cluster.ids <- c(
  "0" = "CD4 T cells",
  "1" = "CD14 Monocytes",
  "2" = "CD4 T cells",
  "3" = "B cells",
  "4" = "CD8 T cells",
  "5" = "NK cells",
  "6" = "CD14 Monocytes",
  "7" = "CD4 T cells",
  "8" = "B cells"
)

# Apply annotations (adjust based on your cluster numbers)
# First, let's create a function to automatically assign based on marker expression
assign_cell_types <- function(seurat_obj) {
  clusters <- levels(Idents(seurat_obj))
  cell_type_assignments <- character(length(clusters))
  names(cell_type_assignments) <- clusters
  
  for (cluster in clusters) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster)
    
    # Calculate mean expression for each marker set
    scores <- sapply(canonical_markers, function(markers) {
      markers_present <- markers[markers %in% rownames(seurat_obj)]
      if (length(markers_present) > 0) {
        mean(colMeans(GetAssayData(seurat_obj, slot = "data")[markers_present, cluster_cells, drop = FALSE]))
      } else {
        0
      }
    })
    
    # Assign cell type with highest score
    cell_type_assignments[cluster] <- names(which.max(scores))
  }
  
  return(cell_type_assignments)
}

# Automatically assign cell types
cell_type_labels <- assign_cell_types(pbmc)
print("Assigned cell types:")
print(cell_type_labels)

# Rename clusters with cell type labels
names(cell_type_labels) <- levels(Idents(pbmc))
pbmc <- RenameIdents(pbmc, cell_type_labels)

# Store cell type in metadata
pbmc$cell_type <- Idents(pbmc)

# ============================================================================
# STEP 12: Final Visualization
# ============================================================================

# UMAP with cell type labels
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5) + 
  ggtitle("PBMC Cell Types (2,700 cells)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p1)

# UMAP without labels for cleaner view
p2 <- DimPlot(pbmc, reduction = "umap", pt.size = 0.5) + 
  ggtitle("PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p2)

# Cell type composition
cell_type_counts <- table(pbmc$cell_type)
print("Cell type composition:")
print(cell_type_counts)
print(prop.table(cell_type_counts))

# Bar plot of cell type composition
cell_type_df <- as.data.frame(cell_type_counts)
colnames(cell_type_df) <- c("Cell_Type", "Count")

ggplot(cell_type_df, aes(x = reorder(Cell_Type, -Count), y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PBMC Cell Type Distribution",
       x = "Cell Type",
       y = "Number of Cells") +
  geom_text(aes(label = Count), vjust = -0.5)

# ============================================================================
# STEP 13: Save Results
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_annotated.rds")

# Save UMAP plot
ggsave("pbmc_umap_cell_types.pdf", plot = p1, width = 10, height = 8)
ggsave("pbmc_umap_cell_types.png", plot = p1, width = 10, height = 8, dpi = 300)

# Export cell type annotations
cell_annotations <- data.frame(
  barcode = colnames(pbmc),
  cell_type = pbmc$cell_type,
  cluster = pbmc$seurat_clusters
)
write.csv(cell_annotations, "cell_type_annotations.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")
cat("Total cells analyzed:", ncol(pbmc), "\n")
cat("Cell types identified:", length(unique(pbmc$cell_type)), "\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
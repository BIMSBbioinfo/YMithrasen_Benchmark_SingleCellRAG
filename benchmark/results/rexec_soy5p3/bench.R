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

# Read barcodes (cell IDs)
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(barcodes) <- "barcode"

# Read genes (gene IDs and symbols)
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

cat("Seurat object created with", ncol(pbmc), "cells and", nrow(pbmc), "genes\n")

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
# STEP 4: Feature Selection (Highly Variable Genes)
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

# Elbow plot to determine dimensionality
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

# Visualize UMAP with cluster labels
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP - Seurat Clusters")

# ============================================================================
# STEP 9: Find Cluster Markers
# ============================================================================

# Find markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers per cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

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

# Visualize canonical markers on UMAP
for (cell_type in names(canonical_markers)) {
  markers <- canonical_markers[[cell_type]]
  available_markers <- markers[markers %in% rownames(pbmc)]
  
  if (length(available_markers) > 0) {
    print(FeaturePlot(pbmc, features = available_markers, ncol = length(available_markers)) + 
            plot_annotation(title = cell_type))
  }
}

# Visualize all canonical markers together
all_markers <- unlist(canonical_markers)
all_markers <- all_markers[all_markers %in% rownames(pbmc)]

DotPlot(pbmc, features = all_markers) + 
  RotatedAxis() + 
  ggtitle("Canonical Marker Expression by Cluster")

# ============================================================================
# STEP 11: Manual Cell Type Assignment
# ============================================================================

# Calculate average expression of markers per cluster
cluster_avg_exp <- AverageExpression(pbmc, features = all_markers, group.by = "seurat_clusters")

# Manual annotation based on marker expression
# You'll need to inspect the plots and adjust these assignments
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

# Apply cell type labels
current_clusters <- levels(Idents(pbmc))
pbmc$cell_type <- as.character(Idents(pbmc))

for (cluster in names(new.cluster.ids)) {
  if (cluster %in% current_clusters) {
    pbmc$cell_type[pbmc$cell_type == cluster] <- new.cluster.ids[cluster]
  }
}

# Set cell type as active identity
Idents(pbmc) <- pbmc$cell_type

# ============================================================================
# STEP 12: Final Visualization
# ============================================================================

# UMAP with cell type labels
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.8, label.size = 5) + 
  ggtitle("PBMC Cell Types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(p1)

# Save the plot
ggsave("PBMC_UMAP_CellTypes.pdf", plot = p1, width = 10, height = 8)
ggsave("PBMC_UMAP_CellTypes.png", plot = p1, width = 10, height = 8, dpi = 300)

# Cell type composition
cell_type_counts <- table(pbmc$cell_type)
print("Cell type composition:")
print(cell_type_counts)

# Proportions
cell_type_props <- prop.table(cell_type_counts) * 100
print("\nCell type proportions (%):")
print(round(cell_type_props, 2))

# Bar plot of cell type composition
p2 <- ggplot(data.frame(cell_type = names(cell_type_counts), 
                        count = as.numeric(cell_type_counts)), 
             aes(x = reorder(cell_type, -count), y = count, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "PBMC Cell Type Composition", x = "Cell Type", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p2)

# ============================================================================
# STEP 13: Marker Heatmap
# ============================================================================

# Heatmap of canonical markers
DoHeatmap(pbmc, features = all_markers, size = 4) + 
  ggtitle("Canonical Marker Expression Heatmap")

# ============================================================================
# STEP 14: Save Results
# ============================================================================

# Save the Seurat object
saveRDS(pbmc, file = "pbmc_seurat_annotated.rds")

# Save cell type annotations
write.csv(data.frame(barcode = colnames(pbmc), 
                     cluster = pbmc$seurat_clusters,
                     cell_type = pbmc$cell_type),
          file = "cell_type_annotations.csv", 
          row.names = FALSE)

cat("\nAnalysis complete! Results saved.\n")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
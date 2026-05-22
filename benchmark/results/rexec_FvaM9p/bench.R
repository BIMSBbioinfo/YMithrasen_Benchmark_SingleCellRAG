setwd("C:/MDC_internship/rolvdesktop/rolvdesktop_electron/benchmark/datasets/filtered_gene_bc_matrices")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("ggplot2", quietly=TRUE)) tryCatch(BiocManager::install("ggplot2", ask=FALSE, quiet=TRUE), error=function(e) install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(ggplot2)))

# Load the PBMC dataset
sc.data <- Read10X(data.dir = ".")

# Create Seurat object
sc <- CreateSeuratObject(counts = sc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Add mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# QC visualization
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Preprocessing: filter, normalize, find variable features, and scale
sc_preprocessed <- sc %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(.))

# Visualize top 10 variable features
top10 <- head(VariableFeatures(sc_preprocessed), 10)
plot1 <- VariableFeaturePlot(sc_preprocessed)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

# PCA
sc_preprocessed <- RunPCA(sc_preprocessed, features = VariableFeatures(object = sc_preprocessed))

# Clustering
set.seed(42)
sc_preprocessed <- FindNeighbors(sc_preprocessed, dims = 1:10)
sc_preprocessed <- FindClusters(sc_preprocessed, resolution = 0.5)

# UMAP
set.seed(42)
sc_preprocessed <- RunUMAP(sc_preprocessed, dims = 1:10)

# Visualize clusters before annotation
DimPlot(sc_preprocessed, reduction = "umap", label = TRUE)

# Visualize canonical marker genes
FeaturePlot(sc_preprocessed, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7"))

VlnPlot(sc_preprocessed, features = c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "GNLY", "NKG7"))

# Find all markers to help with annotation
sc.markers <- FindAllMarkers(sc_preprocessed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- sc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
print(top_markers)

# Assign cell type identities based on canonical markers
# Adjust cluster numbers based on marker expression patterns observed
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(sc_preprocessed)
sc_preprocessed <- RenameIdents(sc_preprocessed, new.cluster.ids)

# Visualize final UMAP with cell type labels
plot_final <- DimPlot(sc_preprocessed, reduction = "umap", label = TRUE, label.size = 4.5) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 10)))
plot_final
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
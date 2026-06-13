setwd("C:/MDC_internship/YMithrasen_Benchmark_SingleCellRAG/benchmark/datasets/filtered_gene_bc_matrices/hg19")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
if (!requireNamespace("Seurat", quietly=TRUE)) tryCatch(BiocManager::install("Seurat", ask=FALSE, quiet=TRUE), error=function(e) install.packages("Seurat", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("dplyr", quietly=TRUE)) tryCatch(BiocManager::install("dplyr", ask=FALSE, quiet=TRUE), error=function(e) install.packages("dplyr", repos="https://cloud.r-project.org", quiet=TRUE))
if (!requireNamespace("patchwork", quietly=TRUE)) tryCatch(BiocManager::install("patchwork", ask=FALSE, quiet=TRUE), error=function(e) install.packages("patchwork", repos="https://cloud.r-project.org", quiet=TRUE))
tryCatch({
  suppressMessages(suppressWarnings({
# Load required libraries
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

# ============================================
# 1. Load 10x Genomics Data
# ============================================

# Load the data from the 3 files (barcodes.tsv, genes.tsv, matrix.mtx)
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "PBMC_doublet_removal",
                           min.cells = 3,
                           min.features = 200)

# ============================================
# 2. Initial QC and Filtering
# ============================================

# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 6000 & 
                              percent.mt < 15)

# ============================================
# 3. Standard Pre-processing
# ============================================

# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                             nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================
# 4. Doublet Detection with DoubletFinder
# ============================================

# Run UMAP and clustering (required for DoubletFinder)
pbmc <- RunUMAP(pbmc, dims = 1:30)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Pre-process for DoubletFinder
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Visualize pK selection
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_classic()

# Select optimal pK (usually the peak)
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Estimate doublet rate (adjust based on your cell recovery)
# Rule of thumb: ~0.8% per 1000 cells recovered
# For 5000 cells: ~4% doublet rate
nExp_poi <- round(0.04 * ncol(pbmc))  # Adjust 0.04 based on your experiment

# Homotypic doublet proportion estimate
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
pbmc <- doubletFinder_v3(pbmc, 
                         PCs = 1:30, 
                         pN = 0.25, 
                         pK = pK_optimal, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE, 
                         sct = FALSE)

# The column name will be something like "DF.classifications_0.25_0.09_123"
# Find the correct column name
doublet_col <- grep("DF.classifications", colnames(pbmc@meta.data), value = TRUE)

# Visualize doublets
DimPlot(pbmc, reduction = "umap", group.by = doublet_col)

# ============================================
# 5. Remove Doublets
# ============================================

# Count doublets
table(pbmc@meta.data[[doublet_col]])

# Remove doublets
pbmc_singlets <- subset(pbmc, subset = get(doublet_col) == "Singlet")

# Check remaining cells
print(paste("Original cells:", ncol(pbmc)))
print(paste("Singlets retained:", ncol(pbmc_singlets)))
print(paste("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets)))

# ============================================
# 6. Re-process Clean Data
# ============================================

# Re-run standard workflow on clean data
pbmc_singlets <- NormalizeData(pbmc_singlets)
pbmc_singlets <- FindVariableFeatures(pbmc_singlets, nfeatures = 2000)
pbmc_singlets <- ScaleData(pbmc_singlets)
pbmc_singlets <- RunPCA(pbmc_singlets)
pbmc_singlets <- RunUMAP(pbmc_singlets, dims = 1:30)
pbmc_singlets <- FindNeighbors(pbmc_singlets, dims = 1:30)
pbmc_singlets <- FindClusters(pbmc_singlets, resolution = 0.5)

# Visualize clean data
DimPlot(pbmc_singlets, reduction = "umap", label = TRUE)

# ============================================
# 7. Save Results
# ============================================

# Save cleaned Seurat object
saveRDS(pbmc_singlets, file = "pbmc_singlets_cleaned.rds")

# Export doublet information
write.csv(pbmc@meta.data, file = "doublet_detection_results.csv")

# Create comparison plot
p1 <- DimPlot(pbmc, reduction = "umap", group.by = doublet_col) + 
      ggtitle("Before Doublet Removal")
p2 <- DimPlot(pbmc_singlets, reduction = "umap", label = TRUE) + 
      ggtitle("After Doublet Removal")

library(patchwork)
p1 | p2
ggsave("doublet_removal_comparison.pdf", width = 14, height = 6)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
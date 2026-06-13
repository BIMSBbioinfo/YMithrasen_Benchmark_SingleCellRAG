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

# Load the three files (barcodes.tsv, genes.tsv, matrix.mtx)
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "PBMC_doublet_removal",
                           min.cells = 3,
                           min.features = 200)

# ============================================
# 2. Quality Control and Preprocessing
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
# 3. Standard Preprocessing for DoubletFinder
# ============================================

# Normalize data
pbmc <- NormalizeData(pbmc)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# Run UMAP (use appropriate number of PCs based on elbow plot)
pbmc <- RunUMAP(pbmc, dims = 1:30)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ============================================
# 4. DoubletFinder Parameter Optimization
# ============================================

# pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Visualize pK optimization
ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_classic()

# Select optimal pK (highest BCmetric)
optimal_pK <- bcmvn %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) %>% 
  as.numeric()

print(paste("Optimal pK:", optimal_pK))

# ============================================
# 5. Homotypic Doublet Proportion Estimate
# ============================================

# Get annotations (cluster assignments)
annotations <- pbmc@meta.data$seurat_clusters

# Calculate homotypic doublet proportion
homotypic.prop <- modelHomotypic(annotations)

# Estimate doublet rate based on cell recovery
# Typical 10x doublet rates:
# ~3,000 cells: 2.3%
# ~5,000 cells: 3.9%
# ~8,000 cells: 6.1%
# ~10,000 cells: 7.6%

nExp_poi <- round(0.075 * nrow(pbmc@meta.data))  # Assuming ~7.5% doublet rate
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

print(paste("Expected doublets:", nExp_poi))
print(paste("Adjusted for homotypic doublets:", nExp_poi.adj))

# ============================================
# 6. Run DoubletFinder
# ============================================

pbmc <- doubletFinder_v3(pbmc, 
                         PCs = 1:30, 
                         pN = 0.25,  # proportion of artificial doublets
                         pK = optimal_pK, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE,
                         sct = FALSE)

# Get the DoubletFinder classification column name
# (it will be something like "DF.classifications_0.25_0.09_913")
df_column <- colnames(pbmc@meta.data)[grep("DF.classifications", colnames(pbmc@meta.data))]

# Rename for convenience
pbmc@meta.data$doublet_classification <- pbmc@meta.data[[df_column]]

# ============================================
# 7. Visualize Doublets
# ============================================

# UMAP colored by doublet classification
DimPlot(pbmc, reduction = "umap", group.by = "doublet_classification", 
        cols = c("Singlet" = "blue", "Doublet" = "red")) +
  ggtitle("Doublet Detection Results")

# Count doublets
table(pbmc@meta.data$doublet_classification)

# Visualize doublet score distribution
pANN_column <- colnames(pbmc@meta.data)[grep("pANN", colnames(pbmc@meta.data))]
pbmc@meta.data$doublet_score <- pbmc@meta.data[[pANN_column]]

ggplot(pbmc@meta.data, aes(x = doublet_score, fill = doublet_classification)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  theme_classic() +
  ggtitle("Doublet Score Distribution")

# ============================================
# 8. Remove Doublets
# ============================================

# Filter out doublets
pbmc_singlets <- subset(pbmc, subset = doublet_classification == "Singlet")

print(paste("Original cells:", ncol(pbmc)))
print(paste("After doublet removal:", ncol(pbmc_singlets)))
print(paste("Removed:", ncol(pbmc) - ncol(pbmc_singlets), "doublets"))

# ============================================
# 9. Save Results
# ============================================

# Save cleaned Seurat object
saveRDS(pbmc_singlets, file = "pbmc_singlets_cleaned.rds")

# Export doublet information
write.csv(pbmc@meta.data, file = "doublet_detection_results.csv")

# ============================================
# 10. Re-process Cleaned Data (Optional)
# ============================================

# After doublet removal, re-run standard workflow
pbmc_singlets <- NormalizeData(pbmc_singlets)
pbmc_singlets <- FindVariableFeatures(pbmc_singlets)
pbmc_singlets <- ScaleData(pbmc_singlets)
pbmc_singlets <- RunPCA(pbmc_singlets)
pbmc_singlets <- RunUMAP(pbmc_singlets, dims = 1:30)
pbmc_singlets <- FindNeighbors(pbmc_singlets, dims = 1:30)
pbmc_singlets <- FindClusters(pbmc_singlets, resolution = 0.5)

# Visualize cleaned data
DimPlot(pbmc_singlets, reduction = "umap", label = TRUE) +
  ggtitle("PBMC Clusters After Doublet Removal")

print("Doublet detection and removal complete!")
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
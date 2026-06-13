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
# Assuming your files are in a directory called "filtered_feature_bc_matrix"
pbmc.data <- Read10X(data.dir = ".")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "PBMC", 
                           min.cells = 3, 
                           min.features = 200)

# ============================================
# 2. Quality Control and Filtering
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
pbmc <- NormalizeData(pbmc)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data
pbmc <- ScaleData(pbmc)

# Run PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Determine dimensionality
ElbowPlot(pbmc, ndims = 50)

# ============================================
# 4. Clustering for DoubletFinder
# ============================================
# Use appropriate number of PCs (adjust based on ElbowPlot)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:30)

# ============================================
# 5. DoubletFinder - Parameter Optimization
# ============================================
# pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select optimal pK
pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
cat("Optimal pK:", pK_optimal, "\n")

# ============================================
# 6. Estimate Doublet Rate
# ============================================
# Homotypic doublet proportion estimate
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# Expected doublet rate (adjust based on cell recovery)
# Rule of thumb: ~0.8% per 1000 cells recovered
nCells <- ncol(pbmc)
expected_doublet_rate <- 0.008 * (nCells / 1000)  # Adjust as needed
cat("Expected doublet rate:", expected_doublet_rate, "\n")

# Adjust for homotypic doublets
nExp_poi <- round(expected_doublet_rate * nCells)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

cat("Expected doublets:", nExp_poi, "\n")
cat("Adjusted expected doublets:", nExp_poi.adj, "\n")

# ============================================
# 7. Run DoubletFinder
# ============================================
pbmc <- doubletFinder_v3(pbmc, 
                         PCs = 1:30, 
                         pN = 0.25,  # proportion of artificial doublets
                         pK = pK_optimal, 
                         nExp = nExp_poi.adj,
                         reuse.pANN = FALSE,
                         sct = FALSE)

# Get the DoubletFinder column name (it's auto-generated)
df_column <- colnames(pbmc@meta.data)[grep("DF.classifications", colnames(pbmc@meta.data))]

# ============================================
# 8. Visualize Doublets
# ============================================
# UMAP with doublet labels
DimPlot(pbmc, reduction = "umap", group.by = df_column, 
        cols = c("Singlet" = "blue", "Doublet" = "red")) +
  ggtitle("Doublet Detection Results")

# Count doublets
table(pbmc@meta.data[[df_column]])

# Doublets per cluster
doublet_summary <- pbmc@meta.data %>%
  group_by(seurat_clusters, .data[[df_column]]) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

print(doublet_summary)

# ============================================
# 9. Remove Doublets
# ============================================
# Filter out doublets
pbmc_singlets <- subset(pbmc, subset = .data[[df_column]] == "Singlet")

cat("\nOriginal cells:", ncol(pbmc), "\n")
cat("Singlets retained:", ncol(pbmc_singlets), "\n")
cat("Doublets removed:", ncol(pbmc) - ncol(pbmc_singlets), "\n")

# ============================================
# 10. Re-process Cleaned Data
# ============================================
# Re-run standard workflow on cleaned data
pbmc_singlets <- NormalizeData(pbmc_singlets)
pbmc_singlets <- FindVariableFeatures(pbmc_singlets)
pbmc_singlets <- ScaleData(pbmc_singlets)
pbmc_singlets <- RunPCA(pbmc_singlets)
pbmc_singlets <- FindNeighbors(pbmc_singlets, dims = 1:30)
pbmc_singlets <- FindClusters(pbmc_singlets, resolution = 0.5)
pbmc_singlets <- RunUMAP(pbmc_singlets, dims = 1:30)

# Visualize cleaned data
DimPlot(pbmc_singlets, reduction = "umap", label = TRUE) +
  ggtitle("PBMC Data After Doublet Removal")

# ============================================
# 11. Save Results
# ============================================
# Save cleaned Seurat object
saveRDS(pbmc_singlets, file = "pbmc_singlets_cleaned.rds")

# Export doublet information
write.csv(pbmc@meta.data, file = "doublet_detection_results.csv")

# Save plots
ggsave("doublet_umap.pdf", width = 8, height = 6)
  }))
  cat("BENCH_EXEC_OK\n")
}, error = function(e) {
  cat("BENCH_EXEC_ERROR:", conditionMessage(e), "\n", file=stderr())
})
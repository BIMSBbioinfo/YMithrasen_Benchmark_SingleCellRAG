library(Seurat)
library(patchwork)

# =============================================================================
# PARAMETERS
# =============================================================================
mt_pattern   <- "^MT-"   # mitochondrial gene pattern: "^MT-" human, "^mt-" mouse
ribo_pattern <- "^RP[SL]" # ribosomal gene pattern (optional)
compute_ribo <- FALSE     # set TRUE to also compute ribosomal percentage
# =============================================================================

# Mitochondrial percentage
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)

# Ribosomal percentage (optional)
if (compute_ribo) {
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = ribo_pattern)
}

cat("QC metrics computed. Summary:\n")
print(summary(seu@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]))

# Violin plot — inspect to choose filtering thresholds
features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
if (compute_ribo) features_to_plot <- c(features_to_plot, "percent.ribo")
VlnPlot(seu, features = features_to_plot, ncol = length(features_to_plot))

# Scatter plots — relationships between QC metrics
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

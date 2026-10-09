library(Seurat)

# PARAMETERS
n_variable       <- 2000   # number of highly variable features to select
selection_method <- "vst"  # "vst" (default), "mean.var.plot", or "dispersion"

# Suggested ranges:
#   n_variable: 1000–5000 (2000 is standard)


seu_filtered <- FindVariableFeatures(
  seu_filtered,
  selection.method = selection_method,
  nfeatures        = n_variable
)

top10 <- head(VariableFeatures(seu_filtered), 10)
cat("Top 10 variable genes:", paste(top10, collapse = ", "), "\n")

# Plot variable features
plot1 <- VariableFeaturePlot(seu_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,
                     xnudge = 0, ynudge = 0)
print(plot2)

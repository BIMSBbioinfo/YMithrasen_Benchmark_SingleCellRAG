library(Seurat)
library(DoubletFinder)
library(ggplot2)

# PARAMETERS
n_pcs        <- 20     # must match the value used in 10_pca.R
doublet_rate <- 0.07   # expected proportion of doublets (0.05–0.10)
pK_fixed     <- 0.01   # fixed pK value — see note below
run_pk_sweep <- FALSE  # set TRUE to run full pK parameter sweep (slower)

# Note on pK: the optimal pK changes every run because it depends on the data.
# The column name DoubletFinder creates also encodes pK, so it changes too.
# pK = 0.01 is a stable default. For publication results, set run_pk_sweep = TRUE,
# inspect the BCmetric plot, and read the peak pK value from the plot.


# Optional pK sweep
if (run_pk_sweep) {
  sweep.res  <- paramSweep(seu_filtered, PCs = 1:n_pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn      <- find.pK(sweep.stats)
  ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() + geom_line() +
    ggtitle("pK selection — pick value at peak BCmetric")
  pK_fixed <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat("Optimal pK selected:", pK_fixed, "\n")
}

# Homotypic doublet proportion
annotations    <- seu_filtered@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(annotations)

# Expected doublet counts
n_exp          <- round(doublet_rate * nrow(seu_filtered@meta.data))
n_exp_adj      <- round(n_exp * (1 - homotypic_prop))

cat("Expected doublets (raw):     ", n_exp, "\n")
cat("Expected doublets (adjusted):", n_exp_adj, "\n")

# Run DoubletFinder
seu_df <- doubletFinder(seu_filtered, PCs = 1:n_pcs, pN = 0.25, pK = pK_fixed,
                        nExp = n_exp_adj, reuse.pANN = FALSE, sct = FALSE)

# Dynamically find the classification column — never hardcode this
df_col <- grep("DF.classifications", colnames(seu_df@meta.data), value = TRUE)[1]
cat("DoubletFinder column:", df_col, "\n")
print(table(seu_df@meta.data[[df_col]]))

# Visualise doublets on UMAP
DimPlot(seu_df, reduction = "umap", group.by = df_col) +
  ggtitle("DoubletFinder classifications")

# Subset to singlets
seu_singlets <- subset(
  seu_df,
  cells = rownames(seu_df@meta.data)[seu_df@meta.data[[df_col]] == "Singlet"]
)
cat("Singlets retained:", ncol(seu_singlets), "\n")
cat("Doublets removed: ", ncol(seu_df) - ncol(seu_singlets), "\n")

DimPlot(seu_singlets, reduction = "umap", label = TRUE) +
  ggtitle("UMAP — singlets only")

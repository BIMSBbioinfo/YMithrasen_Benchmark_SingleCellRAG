library(Seurat)
library(ggplot2)

# PARAMETERS
n_pcs        <- 10     # must match value used in 10_pca.R
doublet_rate <- 0.07   # expected proportion of doublets (0.05-0.10)
pK_fixed     <- 0.01   # stable default; use sweep only if run_pk_sweep = TRUE
run_pk_sweep <- FALSE  # set TRUE to run full pK sweep (slower)

# ============================================================
# REQUIRED INPUT: seu_filtered
# Must already have: NormalizeData -> FindVariableFeatures ->
# ScaleData -> RunPCA -> FindNeighbors -> FindClusters -> RunUMAP
#
# DO NOT call doubletFinder() directly.
# ALWAYS use run_doubletfinder_safe() defined below.
# ============================================================

run_doubletfinder_safe <- function(seu, n_pcs = 10, doublet_rate = 0.07,
                                   pK = 0.01, run_pk_sweep = FALSE) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    message("DoubletFinder not installed. Marking all cells as Singlet.")
    seu$doublet_status <- "Singlet"
    return(seu)
  }

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    message("No seurat_clusters found. Marking all cells as Singlet.")
    seu$doublet_status <- "Singlet"
    return(seu)
  }

  # Optional pK sweep — only when run_pk_sweep = TRUE
  if (isTRUE(run_pk_sweep)) {
    sweep.res   <- DoubletFinder::paramSweep(seu, PCs = 1:n_pcs, sct = FALSE)
    sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
    bcmvn       <- DoubletFinder::find.pK(sweep.stats)
    # CRITICAL: extract a single numeric scalar from the data frame
    # Never pass bcmvn or bcmvn$pK directly to doubletFinder — it will crash
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("Optimal pK: ", pK)
  }
  # If run_pk_sweep = FALSE, use the pK value passed as argument (default 0.01)

  annotations    <- as.character(seu@meta.data$seurat_clusters)
  homotypic_prop <- tryCatch(
    DoubletFinder::modelHomotypic(annotations),
    error = function(e) { message("modelHomotypic failed: ", e$message); 0 }
  )

  n_exp     <- round(doublet_rate * ncol(seu))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))
  message("Expected doublets (raw):      ", n_exp)
  message("Expected doublets (adjusted): ", n_exp_adj)

  seu_df <- tryCatch(
    DoubletFinder::doubletFinder(
      seu,
      PCs        = 1:n_pcs,
      pN         = 0.25,
      pK         = as.numeric(pK),   # must be a single numeric scalar
      nExp       = n_exp_adj,
      reuse.pANN = NULL,
      sct        = FALSE
    ),
    error = function(e) {
      message("doubletFinder() failed: ", e$message)
      message("Fallback: marking all cells as Singlet.")
      seu$doublet_status <- "Singlet"
      return(seu)
    }
  )

  # Dynamic column detection — never hardcode the column name
  df_col <- grep("DF.classifications", colnames(seu_df@meta.data), value = TRUE)
  df_col <- if (length(df_col) > 0) df_col[1] else NA

  if (is.na(df_col)) {
    message("No DF.classifications column found. Using fallback.")
    seu_df$doublet_status <- "Singlet"
  } else {
    seu_df$doublet_status <- as.character(seu_df@meta.data[[df_col]])
    print(table(seu_df$doublet_status))
  }

  return(seu_df)
}

# ============================================================
# USAGE — always call run_doubletfinder_safe(), never doubletFinder() directly
# pK sweep is disabled by default (run_pk_sweep = FALSE) for speed and stability
# ============================================================
seu_df <- run_doubletfinder_safe(
  seu          = seu_filtered,
  n_pcs        = n_pcs,
  doublet_rate = doublet_rate,
  pK           = pK_fixed,
  run_pk_sweep = run_pk_sweep
)

if ("umap" %in% names(seu_df@reductions)) {
  print(DimPlot(seu_df, reduction = "umap", group.by = "doublet_status") +
          ggtitle("Doublet status"))
}

seu_singlets <- subset(
  seu_df,
  cells = colnames(seu_df)[seu_df$doublet_status == "Singlet"]
)
cat("Singlets retained:", ncol(seu_singlets), "\n")
cat("Doublets removed: ", ncol(seu_df) - ncol(seu_singlets), "\n")

saveRDS(seu_df,       "seurat_doublet_annotated.rds")
saveRDS(seu_singlets, "seurat_singlets.rds")
write.csv(
  data.frame(cell = colnames(seu_df), doublet_status = seu_df$doublet_status),
  "doublet_calls.csv", row.names = FALSE
)
cat("Saved: seurat_doublet_annotated.rds\n")
cat("Saved: seurat_singlets.rds\n")
cat("Saved: doublet_calls.csv\n")
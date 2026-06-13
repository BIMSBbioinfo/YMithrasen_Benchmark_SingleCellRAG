library(Seurat)

# =============================================================================
# PARAMETERS
# =============================================================================
scale_all_genes    <- FALSE          # TRUE scales all genes; FALSE scales only variable features
vars_to_regress    <- NULL           # covariates to regress out, e.g. c("percent.mt", "nCount_RNA")
                                     # set NULL to skip regression
# Examples:
#   vars_to_regress <- c("percent.mt")               # regress mitochondrial content
#   vars_to_regress <- c("S.Score", "G2M.Score")     # regress cell cycle (if not done in 06)
#   vars_to_regress <- c("percent.mt", "nCount_RNA")  # regress multiple covariates
# =============================================================================

features_to_scale <- if (scale_all_genes) rownames(seu_filtered) else VariableFeatures(seu_filtered)

seu_filtered <- ScaleData(
  seu_filtered,
  features        = features_to_scale,
  vars.to.regress = vars_to_regress
)

cat("Scaling complete\n")
cat("Features scaled:", length(features_to_scale), "\n")
if (!is.null(vars_to_regress)) {
  cat("Covariates regressed:", paste(vars_to_regress, collapse = ", "), "\n")
}

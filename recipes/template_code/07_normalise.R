library(Seurat)

# PARAMETERS
norm_method  <- "lognorm"  # "lognorm" for log-normalisation, "sctransform" for SCTransform
scale_factor <- 10000      # scale factor for log-normalisation (standard: 10000)


if (norm_method == "lognorm") {
  seu_filtered <- NormalizeData(
    seu_filtered,
    normalization.method = "LogNormalize",
    scale.factor         = scale_factor
  )
  cat("Log-normalisation complete\n")

} else if (norm_method == "sctransform") {
  # SCTransform: regularised negative binomial regression
  # Removes variation due to sequencing depth without a separate ScaleData step
  seu_filtered <- SCTransform(seu_filtered, verbose = FALSE)
  cat("SCTransform normalisation complete\n")
  cat("Note: ScaleData and FindVariableFeatures are incorporated in SCTransform\n")
  cat("      Skip scripts 08 and 09 if using SCTransform\n")

} else {
  stop("norm_method must be 'lognorm' or 'sctransform'")
}

library(Seurat)
library(dplyr)

# PARAMETERS
integration_method <- "cca"   # "none" (merge only), "cca" (Seurat CCA/RPCA), "harmony"
n_pcs              <- 20      # PCs to use for integration and downstream steps
resolution         <- 0.5     # clustering resolution post-integration
reduction_type     <- "rpca"  # for CCA: "cca" or "rpca" (RPCA is faster for large datasets)
min_common_genes   <- 5000    # warn if fewer than this many genes are shared across samples

# Harmony-specific (only used if integration_method == "harmony")
harmony_group_var  <- "orig.ident"  # metadata column to correct by (e.g. sample, batch, donor)


#  Step 1: Find common gene set across all samples
# Samples from different protocols, references, or library prep batches may
# not share all genes. Integration requires a common feature set.

cat("── Gene set compatibility check ──\n")
gene_counts <- sapply(sample_list, nrow)
cat("Genes per sample:\n")
print(gene_counts)

common_genes <- Reduce(intersect, lapply(sample_list, rownames))
cat("Common genes across all samples:", length(common_genes), "\n")

n_lost <- max(gene_counts) - length(common_genes)
if (n_lost > 0) {
  warning(paste0(
    n_lost, " genes removed because they are absent in one or more samples.\n",
    "  This is normal when samples come from different Cell Ranger references or species.\n",
    "  Check that all samples use the same genome reference if this number is large."
  ))
}

if (length(common_genes) < min_common_genes) {
  warning(paste0(
    "Only ", length(common_genes), " common genes found across samples.\n",
    "  This is fewer than the recommended minimum of ", min_common_genes, ".\n",
    "  Check that samples are from the same species and were aligned to the same reference."
  ))
}

# Subset all samples to common genes
sample_list <- lapply(sample_list, function(seu) {
  subset(seu, features = common_genes)
})
cat("All samples subset to", length(common_genes), "common genes.\n\n")

# Step 2: Preprocess each sample individually
preprocess_sample <- function(seu) {
  seu %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
}

cat("Preprocessing each sample individually...\n")
sample_list_processed <- lapply(sample_list, preprocess_sample)
cat("Preprocessing complete.\n\n")

# Step 3: Integrate
cat("Integration method:", integration_method, "\n")

if (integration_method == "none") {
  cat("Merging without batch correction...\n")
  merged <- merge(
    x          = sample_list_processed[[1]],
    y          = sample_list_processed[-1],
    add.cell.ids = names(sample_list_processed)
  )
  integrated <- merged %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:n_pcs, verbose = FALSE) %>%
    FindNeighbors(dims = 1:n_pcs, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)

  cat("Note: no batch correction applied — samples may cluster by batch.\n")
  cat("      Use group.by = 'orig.ident' in DimPlot to check batch separation.\n")

} else if (integration_method == "cca") {
  cat("Finding integration anchors (", reduction_type, ")...\n")
  anchors <- FindIntegrationAnchors(
    object.list = sample_list_processed,
    dims        = 1:n_pcs,
    reduction   = reduction_type
  )
  cat("Integrating data...\n")
  integrated <- IntegrateData(anchorset = anchors, dims = 1:n_pcs) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:n_pcs, verbose = FALSE) %>%
    FindNeighbors(dims = 1:n_pcs, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)

  DefaultAssay(integrated) <- "integrated"
  cat("Default assay set to 'integrated'.\n")
  cat("Switch to 'RNA' assay for differential expression and marker finding.\n")

} else if (integration_method == "harmony") {
  library(harmony)
  cat("Merging samples for Harmony integration...\n")
  merged <- merge(
    x          = sample_list_processed[[1]],
    y          = sample_list_processed[-1],
    add.cell.ids = names(sample_list_processed)
  ) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE)

  cat("Running Harmony batch correction on variable:", harmony_group_var, "\n")
  integrated <- RunHarmony(merged, group.by.vars = harmony_group_var, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:n_pcs, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", dims = 1:n_pcs, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)

  cat("Harmony reduction stored as 'harmony'.\n")
  cat("Use reduction = 'harmony' for downstream steps.\n")

} else {
  stop("integration_method must be 'none', 'cca', or 'harmony'")
}

# Step 4: Report and visualise
cat("\n── Integration complete ──\n")
cat("Method:   ", integration_method, "\n")
cat("Cells:    ", ncol(integrated), "\n")
cat("Clusters: ", length(levels(Idents(integrated))), "\n")

DimPlot(integrated, reduction = "umap", group.by = "orig.ident") +
  ggtitle(paste0("UMAP by sample — ", integration_method))

DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP by cluster")

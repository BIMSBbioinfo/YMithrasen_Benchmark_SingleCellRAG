library(Seurat)

# PARAMETERS
annotation_method <- "manual"   # "manual" or "singler"

# Manual annotation
# Map cluster numbers (as strings) to cell type labels
# Replace with labels appropriate for your dataset
manual_labels <- c(
  "0" = "Cell type A",
  "1" = "Cell type B",
  "2" = "Cell type C"
  # add one entry per cluster
)

# SingleR automated annotation
singler_species  <- "human"   # "human" or "mouse"
singler_ref      <- "HPCA"    # reference dataset: "HPCA" (human), "ImmGen" (mouse),
                               # or "BlueprintEncode", "MouseRNAseqData", etc.


if (annotation_method == "manual") {
  names(manual_labels) <- levels(seu_filtered)
  seu_filtered <- RenameIdents(seu_filtered, manual_labels)

  DimPlot(seu_filtered, reduction = "umap", label = TRUE, label.size = 4.5) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("UMAP — manual cell type annotation") +
    theme(axis.title       = element_text(size = 18),
          legend.text      = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 10)))

} else if (annotation_method == "singler") {
  library(SingleR)
  library(celldex)

  # Load reference dataset
  if (singler_ref == "HPCA") {
    ref <- HumanPrimaryCellAtlasData()
  } else if (singler_ref == "BlueprintEncode") {
    ref <- BlueprintEncodeData()
  } else if (singler_ref == "ImmGen") {
    ref <- ImmGenData()
  } else if (singler_ref == "MouseRNAseqData") {
    ref <- MouseRNAseqData()
  } else {
    stop("singler_ref not recognised. Use: HPCA, BlueprintEncode, ImmGen, MouseRNAseqData")
  }

  # Join layers if using Seurat v5
  if (inherits(seu_filtered[["RNA"]], "Assay5")) {
    seu_filtered[["RNA"]] <- JoinLayers(seu_filtered[["RNA"]])
  }

  DefaultAssay(seu_filtered) <- "RNA"
  rna_matrix <- GetAssayData(seu_filtered, layer = "data")

  singler_results <- SingleR(
    test   = rna_matrix,
    ref    = ref,
    labels = ref$label.main
  )

  seu_filtered <- AddMetaData(seu_filtered,
                              metadata = singler_results$labels,
                              col.name = "SingleR_Labels")

  DimPlot(seu_filtered, reduction = "umap", group.by = "SingleR_Labels",
          label = TRUE) +
    ggtitle("UMAP — SingleR cell type annotation")

} else {
  stop("annotation_method must be 'manual' or 'singler'")
}

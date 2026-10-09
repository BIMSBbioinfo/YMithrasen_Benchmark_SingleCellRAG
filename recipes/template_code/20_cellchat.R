library(CellChat)
library(Seurat)

# PARAMETERS
cell_label_col <- "manual_labels"   # metadata column with cell type labels
condition_col  <- "orig.ident"      # metadata column to split by condition
conditions     <- NULL              # conditions to compare, e.g. c("BMMC", "CD34")
                                    # set NULL to analyse all cells together
species        <- "human"           # "human" or "mouse"
min_cells_per_type <- 10           # warn if a cell type has fewer cells than this


# Step 1: Assay format check and conversion
# CellChat requires Seurat v3-format Assay objects.
# Seurat v5 introduced Assay5 which is not compatible.
# We detect and convert automatically, and report what was done.

cat("── Assay format check ──\n")
cat("Current assay:", DefaultAssay(integrated), "\n")
cat("Assay class:  ", class(integrated[[DefaultAssay(integrated)]]), "\n")

if (inherits(integrated[["RNA"]], "Assay5")) {
  cat("Seurat v5 Assay5 detected. Converting RNA assay to v3 format for CellChat...\n")
  integrated[["RNA_v3"]] <- as(object = integrated[["RNA"]], Class = "Assay")
  DefaultAssay(integrated) <- "RNA_v3"
  assay_use <- "RNA_v3"
  cat("Conversion complete. Using assay: RNA_v3\n")
  cat("Note: RNA_v3 is a copy of RNA in Seurat v3 format — original RNA assay is unchanged.\n")
} else {
  DefaultAssay(integrated) <- "RNA"
  assay_use <- "RNA"
  cat("Seurat v3 Assay detected. No conversion needed. Using assay: RNA\n")
}

# Step 2: Validate cell type labels
cat("\n── Cell type label check ──\n")
if (!cell_label_col %in% colnames(integrated@meta.data)) {
  stop(paste0(
    "Column '", cell_label_col, "' not found in metadata.\n",
    "  Available columns: ", paste(colnames(integrated@meta.data), collapse = ", "), "\n",
    "  Set cell_label_col to the correct column name."
  ))
}

cell_type_counts <- table(integrated@meta.data[[cell_label_col]])
cat("Cell type distribution:\n")
print(cell_type_counts)

small_types <- names(cell_type_counts[cell_type_counts < min_cells_per_type])
if (length(small_types) > 0) {
  warning(paste0(
    "The following cell types have fewer than ", min_cells_per_type, " cells and ",
    "may produce unreliable communication estimates:\n  ",
    paste(small_types, collapse = ", "), "\n",
    "  Consider merging rare types or increasing min_cells_per_type."
  ))
}

# Step 3: Load CellChat database
cat("\n── Loading CellChat database ──\n")
db <- if (species == "human") {
  cat("Using human CellChat database\n")
  CellChatDB.human
} else if (species == "mouse") {
  cat("Using mouse CellChat database\n")
  CellChatDB.mouse
} else {
  stop("species must be 'human' or 'mouse'")
}

# Step 4: Helper — create and run CellChat object
run_cellchat <- function(seu_sub, db, label_col, label = "") {
  cat("Building CellChat object", if (label != "") paste0("(", label, ")"), "...\n")

  data_mat  <- GetAssayData(seu_sub, assay = assay_use, slot = "counts")
  meta_data <- seu_sub@meta.data

  # Align barcodes
  common_bc <- intersect(colnames(data_mat), rownames(meta_data))
  n_dropped  <- ncol(data_mat) - length(common_bc)
  if (n_dropped > 0) {
    warning(paste(n_dropped, "cells dropped due to barcode mismatch between assay and metadata."))
  }
  data_mat  <- data_mat[, common_bc]
  meta_data <- meta_data[common_bc, , drop = FALSE]

  cat("  Cells:", ncol(data_mat), "| Genes:", nrow(data_mat), "\n")

  cc <- createCellChat(object = data_mat, meta = meta_data, group.by = label_col)
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, population.size = FALSE)
  cc <- filterCommunication(cc)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  cat("  Pathways identified:", length(cc@netP$pathways), "\n")
  return(cc)
}

# Step 5: Run CellChat
if (is.null(conditions)) {
  cat("\n── Running CellChat on all cells ──\n")
  cc_obj <- run_cellchat(integrated, db, cell_label_col)

  netVisual_heatmap(cc_obj, measure = "count",  color.heatmap = "Reds",
                    title.name = "Number of interactions")
  netVisual_heatmap(cc_obj, measure = "weight", color.heatmap = "Reds",
                    title.name = "Interaction strength")

} else {
  # Validate conditions exist
  available <- unique(integrated@meta.data[[condition_col]])
  missing   <- setdiff(conditions, available)
  if (length(missing) > 0) {
    stop(paste0(
      "Conditions not found in '", condition_col, "': ",
      paste(missing, collapse = ", "), "\n",
      "  Available values: ", paste(available, collapse = ", ")
    ))
  }

  cat("\n── Running CellChat per condition ──\n")
  cc_list <- lapply(conditions, function(cond) {
    seu_sub <- subset(integrated,
                      cells = colnames(integrated)[integrated@meta.data[[condition_col]] == cond])
    run_cellchat(seu_sub, db, cell_label_col, label = cond)
  })
  names(cc_list) <- conditions

  # Common pathways
  pathway_lists   <- lapply(cc_list, function(cc) cc@netP$pathways)
  common_pathways <- Reduce(intersect, pathway_lists)
  cat("\nPathways per condition:\n")
  for (cond in conditions) cat(" ", cond, ":", length(pathway_lists[[cond]]), "\n")
  cat("Common pathways:", length(common_pathways), "\n")
  if (length(common_pathways) > 0) {
    cat(" ", paste(common_pathways, collapse = ", "), "\n")
  }

  # Heatmaps per condition
  for (cond in conditions) {
    netVisual_heatmap(cc_list[[cond]], measure = "count",
                      color.heatmap = "Reds",
                      title.name = paste("Interactions —", cond))
    netVisual_heatmap(cc_list[[cond]], measure = "weight",
                      color.heatmap = "Reds",
                      title.name = paste("Interaction strength —", cond))
  }

  # Visualise first common pathway
  if (length(common_pathways) > 0) {
    selected_pathway <- common_pathways[1]
    cat("\nVisualising pathway:", selected_pathway, "\n")
    for (cond in conditions) {
      netVisual_aggregate(cc_list[[cond]], signaling = selected_pathway, layout = "circle")
      title(paste(selected_pathway, "— circle —", cond))
      netVisual_aggregate(cc_list[[cond]], signaling = selected_pathway, layout = "chord")
      title(paste(selected_pathway, "— chord —", cond))
    }
  } else {
    cat("No common pathways found across conditions. Check cell type label consistency.\n")
  }
}

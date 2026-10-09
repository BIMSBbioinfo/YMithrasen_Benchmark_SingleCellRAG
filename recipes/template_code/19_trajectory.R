# PARAMETERS
assay_for_trajectory <- "integrated"  # assay to use: "integrated" or "RNA"
use_partition        <- FALSE          # FALSE = single trajectory across all cells
root_method          <- "auto"         # "auto" = most-connected node;
                                       # "manual" = set root_cell below
root_cell            <- NULL           # cell barcode to use as root when root_method = "manual"
min_cells_trajectory <- 200           # warn if fewer cells than this — trajectory may be unreliable


# Step 1: Dependency check
# Monocle3 and SeuratWrappers must both be installed.
# SeuratWrappers provides as.cell_data_set() for Seurat → Monocle3 conversion.

cat("── Dependency check ──\n")
for (pkg in c("monocle3", "SeuratWrappers")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0(
      "Package '", pkg, "' is not installed.\n",
      "  Install with:\n",
      if (pkg == "monocle3")
        "    remotes::install_github('cole-trapnell-lab/monocle3')\n"
      else
        "    remotes::install_github('satijalab/seurat-wrappers')\n"
    ))
  }
  cat(pkg, ": installed\n")
}

library(monocle3)
library(SeuratWrappers)
library(Seurat)

# Step 2: Cell count check
cat("\n── Dataset check ──\n")
cat("Total cells:", ncol(integrated), "\n")
cat("Assay for trajectory:", assay_for_trajectory, "\n")

if (!assay_for_trajectory %in% names(integrated@assays)) {
  stop(paste0(
    "Assay '", assay_for_trajectory, "' not found in Seurat object.\n",
    "  Available assays: ", paste(names(integrated@assays), collapse = ", "), "\n",
    "  Set assay_for_trajectory to one of the above."
  ))
}

if (ncol(integrated) < min_cells_trajectory) {
  warning(paste0(
    "Only ", ncol(integrated), " cells found. Trajectory inference with fewer than ",
    min_cells_trajectory, " cells may produce unreliable results.\n",
    "  Consider using a larger dataset or interpreting results with caution."
  ))
}

# Step 3: Convert Seurat → Monocle3 CDS
cat("\n── Converting Seurat object to Monocle3 CDS ──\n")
cat("This uses as.cell_data_set() from SeuratWrappers.\n")
cat("UMAP coordinates from Seurat will be transferred to the CDS.\n")

cds <- as.cell_data_set(integrated, assay = assay_for_trajectory)
fData(cds)$gene_short_name <- rownames(fData(cds))
cat("CDS created:", ncol(cds), "cells,", nrow(cds), "genes\n")

# Step 4: Set partition
cat("\n── Setting partition ──\n")
if (!use_partition) {
  cat("use_partition = FALSE: treating all cells as a single trajectory.\n")
  cat("Set use_partition = TRUE if you expect multiple disconnected lineages.\n")
  partition_vec <- rep(1, ncol(cds))
  names(partition_vec) <- colnames(cds)
  cds@clusters$UMAP$partitions <- as.factor(partition_vec)
} else {
  cat("use_partition = TRUE: using Monocle3 partitions.\n")
}

# Step 5: Learn trajectory graph
cat("\n── Learning trajectory graph ──\n")
cds <- learn_graph(cds, use_partition = use_partition, verbose = FALSE)
cat("Trajectory graph learned.\n")

# Step 6: Order cells by pseudotime
cat("\n── Ordering cells by pseudotime ──\n")

if (root_method == "auto") {
  cat("root_method = 'auto': selecting root as the most-connected principal graph node.\n")
  cat("This selects the node with the highest number of neighbouring cells as the start point.\n")
  cat("If results look biologically wrong, switch to root_method = 'manual' and specify root_cell.\n")

  get_earliest_node <- function(cds) {
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    node_freq      <- table(closest_vertex)
    root_node      <- paste0("Y_", names(which.max(node_freq)))
    return(root_node)
  }

  root_node <- get_earliest_node(cds)
  cat("Root node selected:", root_node, "\n")

} else if (root_method == "manual") {
  if (is.null(root_cell)) {
    stop(paste0(
      "root_method = 'manual' requires root_cell to be set.\n",
      "  Set root_cell to a valid cell barcode from your dataset.\n",
      "  Example: root_cell <- colnames(integrated)[1]"
    ))
  }
  if (!root_cell %in% colnames(cds)) {
    stop(paste0(
      "root_cell '", root_cell, "' not found in the CDS.\n",
      "  Check that it is a valid cell barcode."
    ))
  }
  root_node <- root_cell
  cat("Root cell set manually:", root_node, "\n")

} else {
  stop("root_method must be 'auto' or 'manual'")
}

cds <- order_cells(cds, root_pr_nodes = root_node)
cat("Pseudotime computed.\n")

# Step 7: Visualise
cat("\n── Plotting trajectory ──\n")
plot_cells(
  cds,
  color_cells_by          = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_branch_points     = FALSE,
  label_roots             = FALSE,
  label_leaves            = FALSE
) + ggtitle("Trajectory — pseudotime")

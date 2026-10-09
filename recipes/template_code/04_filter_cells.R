library(Seurat)

# PARAMETERS — set based on QC plots from 03_qc_metrics.R
min_feat      <- 200    # minimum genes per cell
max_feat      <- 2500   # maximum genes per cell
min_counts    <- 800    # minimum total UMI counts per cell (set 0 to skip)
max_counts    <- Inf    # maximum total UMI counts (set Inf to skip)
max_mt        <- 5      # maximum mitochondrial percentage
filter_method <- "multiplets"  # "none" or "multiplets" (if multiplets column exists)
multiplet_col <- "multiplets"  # metadata column name for doublet/singlet labels
singlet_label <- "singlet"     # value in multiplet_col that indicates a singlet

# Suggested ranges:
#   min_feat:   100–500  (lower tail of nFeature_RNA violin)
#   max_feat:   2000–6000 (upper tail; complex tissues may need higher)
#   max_mt:     5–25     (inspect violin; solid tissues often 10–25)
#   min_counts: 500–1000


n_before <- ncol(seu)

# Build filter expression
keep <- seu$nFeature_RNA > min_feat &
        seu$nFeature_RNA < max_feat &
        seu$percent.mt   < max_mt

if (min_counts > 0)    keep <- keep & seu$nCount_RNA > min_counts
if (is.finite(max_counts)) keep <- keep & seu$nCount_RNA < max_counts

# Optional: filter pre-annotated multiplets (e.g. from 10X CellRanger or ExperimentHub)
if (filter_method == "multiplets" && multiplet_col %in% colnames(seu@meta.data)) {
  keep <- keep & seu@meta.data[[multiplet_col]] == singlet_label
}

seu_filtered <- subset(seu, cells = colnames(seu)[keep])

cat("Cells before filtering:", n_before, "\n")
cat("Cells after filtering: ", ncol(seu_filtered), "\n")
cat("Cells removed:         ", n_before - ncol(seu_filtered), "\n")

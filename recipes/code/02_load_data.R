library(Seurat)

# =============================================================================
# PARAMETERS
# =============================================================================
input_type   <- "mtx"          # one of: "mtx", "h5", "rds", "sce", "matrix"
data_path    <- "path/to/data" # for mtx: folder path; for others: file path
project_name <- "my_project"   # label for the Seurat object
min_cells    <- 3              # include genes detected in at least this many cells
min_features <- 200            # include cells with at least this many genes

# Suggested ranges:
#   min_cells:    2–5
#   min_features: 100–500
# =============================================================================

if (input_type == "mtx") {
  # 10X Cell Ranger MTX output — handles both v2 (genes.tsv) and v3 (features.tsv)
  counts <- Read10X(data.dir = data_path)
  seu <- CreateSeuratObject(counts = counts, project = project_name,
                            min.cells = min_cells, min.features = min_features)

} else if (input_type == "h5") {
  # 10X Cell Ranger HDF5 output
  counts <- Read10X_h5(filename = data_path)
  seu <- CreateSeuratObject(counts = counts, project = project_name,
                            min.cells = min_cells, min.features = min_features)

} else if (input_type == "rds") {
  # Pre-existing Seurat object saved as RDS
  seu <- readRDS(data_path)
  cat("Loaded existing Seurat object from:", data_path, "\n")

} else if (input_type == "sce") {
  # SingleCellExperiment object — requires SeuratWrappers or as.Seurat
  library(SingleCellExperiment)
  sce <- readRDS(data_path)
  seu <- as.Seurat(sce, data = NULL)

} else if (input_type == "matrix") {
  # Tab-delimited expression matrix (rows = genes, columns = cells)
  exp_mat <- read.table(file = data_path, header = TRUE,
                        as.is = TRUE, row.names = 1)
  seu <- CreateSeuratObject(
    counts       = Matrix::Matrix(as.matrix(exp_mat), sparse = TRUE),
    project      = project_name,
    min.cells    = min_cells,
    min.features = min_features
  )
} else {
  stop("input_type must be one of: mtx, h5, rds, sce, matrix")
}

cat("Seurat object created\n")
cat("Genes: ", nrow(seu), "\n")
cat("Cells: ", ncol(seu), "\n")
print(seu)

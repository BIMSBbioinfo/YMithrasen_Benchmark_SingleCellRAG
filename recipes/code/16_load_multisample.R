library(Seurat)
library(dplyr)

# PARAMETERS
input_type <- "rds"   # "rds", "mtx", or "h5"
data_dir   <- "path/to/samples"  # folder containing all sample files

# Sample metadata table — one row per sample
# Replace with your own sample information
sample_metadata <- data.frame(
  file      = c("sample1.rds", "sample2.rds"),   # file names inside data_dir
  name      = c("Sample1",     "Sample2"),         # sample labels
  condition = c("control",     "treated"),         # experimental condition
  donor     = c("D1",          "D2"),              # donor / subject ID
  replicate = c("T1",          "T1"),              # replicate
  stringsAsFactors = FALSE
  # add more columns as needed: sex, batch, timepoint, etc.
)

min_cells    <- 3
min_features <- 200
mt_pattern   <- "^MT-"   # "^MT-" human, "^mt-" mouse


# Load each sample, create Seurat object, add metadata, compute percent.mt
load_sample <- function(i, sample_metadata, data_dir, input_type,
                        min_cells, min_features) {
  filepath <- file.path(data_dir, sample_metadata$file[i])

  if (input_type == "rds") {
    counts <- readRDS(filepath)
    if (inherits(counts, "Seurat")) return(counts)  # already a Seurat object
  } else if (input_type == "mtx") {
    counts <- Read10X(data.dir = filepath)
  } else if (input_type == "h5") {
    counts <- Read10X_h5(filename = filepath)
  }

  seu <- CreateSeuratObject(counts = counts, project = sample_metadata$name[i],
                            min.cells = min_cells, min.features = min_features)
  return(seu)
}

sample_list <- lapply(seq_len(nrow(sample_metadata)), function(i) {
  seu <- load_sample(i, sample_metadata, data_dir, input_type,
                     min_cells, min_features)

  # Add metadata columns
  for (col in setdiff(colnames(sample_metadata), "file")) {
    seu[[col]] <- sample_metadata[[col]][i]
  }

  # Compute mitochondrial percentage
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)

  cat("Loaded:", sample_metadata$name[i],
      "| Cells:", ncol(seu),
      "| Genes:", nrow(seu), "\n")
  return(seu)
})

names(sample_list) <- sample_metadata$name

# QC violin plots per sample
for (nm in names(sample_list)) {
  print(VlnPlot(sample_list[[nm]],
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) +
          ggtitle(nm))
}

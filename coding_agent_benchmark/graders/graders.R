# graders.R
# Deterministic graders for coding agent benchmark
# Prompt IDs: simple_01..10, medium_01..10, hard_01..10
#
# DESIGN PHILOSOPHY (generalized from the original PBMC-exact-value version):
#   Prefer artifact/state checks (does the resulting Seurat object actually
#   contain the expected structure — a reduction, a metadata column, a
#   non-trivial cluster assignment) over exact dataset-specific numbers
#   (an exact cell count, an exact cluster-count range, an exact row count).
#   Source-code (grepl) checks are used only as a fallback for partial
#   credit when no valid artifact exists to inspect.
#
#   A few checks retain specific values because they plausibly reflect the
#   TASK PROMPT's own instructions rather than an assumption about this
#   particular dataset — e.g. "compare cluster 0 vs cluster 3", "save as
#   pbmc_final.rds", named marker genes the prompt explicitly asks for).
#   These are marked "TASK-SPECIFIED, not dataset-specific — verify against
#   prompts file" below; if the actual prompt text doesn't require them,
#   loosen further.
#
# Called as: Rscript graders.R <prompt_id> <work_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

# ── Helpers ──────────────────────────────────────────────────────────────────

find_rds <- function(work_dir) list.files(work_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

find_r_files <- function(work_dir) {
  files <- list.files(work_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  normalized <- gsub("\\\\", "/", files)
  files[!grepl("/recipes/", normalized, fixed = TRUE)]
}

find_csv <- function(work_dir) list.files(work_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
find_pdf <- function(work_dir) list.files(work_dir, pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)

r_content <- function(work_dir) {
  files <- find_r_files(work_dir)
  if (length(files) == 0) return("")
  paste(sapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = " ")), collapse = " ")
}

load_seurat <- function(work_dir) {
  rds <- find_rds(work_dir)
  for (f in rds) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(obj) && inherits(obj, "Seurat")) return(obj)
  }
  NULL
}

pass    <- function(reason, score = 1.0) list(pass = TRUE,  score = score, reason = reason)
partial <- function(reason, score)       list(pass = FALSE, score = score, reason = reason)
fail    <- function(reason)              list(pass = FALSE, score = 0.0,  reason = reason)

# ── Grounding rate (unchanged — shared with skills_benchmark_runner.js /
#    recipe_benchmark_runner.js, which still consume this field) ──────────────

RECIPE_FUNCTIONS <- c(
  "Read10X", "CreateSeuratObject", "PercentageFeatureSet", "NormalizeData",
  "FindVariableFeatures", "ScaleData", "RunPCA", "ElbowPlot", "FindNeighbors",
  "FindClusters", "RunUMAP", "FindMarkers", "FindAllMarkers", "DoHeatmap",
  "FeaturePlot", "VlnPlot", "CellCycleScoring", "SCTransform", "doubletFinder",
  "SingleR", "RenameIdents", "VariableFeaturePlot", "DimPlot", "FeatureScatter"
)

check_grounding <- function(work_dir) {
  tryCatch({
    content <- r_content(work_dir)
    if (nchar(content) == 0) return(list(rate = 0.0, found = character(0)))
    found <- RECIPE_FUNCTIONS[sapply(RECIPE_FUNCTIONS, function(fn) grepl(fn, content, fixed = TRUE))]
    list(rate = length(found) / length(RECIPE_FUNCTIONS), found = found)
  }, error = function(e) list(rate = 0.0, found = character(0)))
}

# ── Graders ──────────────────────────────────────────────────────────────────

grade_simple_01 <- function(work_dir) {
  # Artifact-first: any valid, non-empty Seurat object. No dataset-specific
  # cell/gene count assumed — a differently-sized dataset would still pass.
  seu <- load_seurat(work_dir)
  if (!is.null(seu) && inherits(seu, "Seurat") && ncol(seu) > 0 && nrow(seu) > 0) {
    return(pass(sprintf("Valid Seurat object created (%d cells, %d features)", ncol(seu), nrow(seu))))
  }
  content <- r_content(work_dir)
  if (grepl("Read10X|Read10X_h5|readMM", content) && grepl("CreateSeuratObject", content)) {
    return(partial("Valid loading workflow found in code but no saved object", 0.7))
  }
  fail("No valid single-cell object creation found")
}

grade_simple_02 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)

  has_mt <- (!is.null(seu) && any(grepl("percent\\.?mt|mito", colnames(seu@meta.data), ignore.case = TRUE))) ||
            grepl("PercentageFeatureSet|percent\\.mt", content)
  has_vlnplot_code <- grepl("VlnPlot", content)
  has_plot_file <- length(find_pdf(work_dir)) > 0

  if (has_mt && has_vlnplot_code && has_plot_file) return(pass("QC metrics and violin plot generated"))
  if (has_mt && has_vlnplot_code) return(partial("VlnPlot called but no output file confirmed", 0.8))
  if (has_mt) return(partial("QC metrics present but no visualisation found", 0.5))
  fail("QC metrics not found")
}

grade_simple_03 <- function(work_dir) {
  # No exact dataset-specific filtered-cell-count range — checks that a
  # valid object exists AND the code shows evidence of an actual filtering
  # step, rather than assuming what post-filter count this dataset "should"
  # produce.
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)
  has_filter_code <- grepl("subset\\s*\\(", content) &&
                      grepl("nFeature_RNA|percent\\.mt|nCount_RNA", content)

  if (!is.null(seu) && ncol(seu) > 0 && has_filter_code) {
    return(pass(sprintf("Filtered Seurat object created (%d cells) with filtering code present", ncol(seu))))
  }
  if (!is.null(seu) && ncol(seu) > 0) {
    return(partial("Valid filtered object present but filtering logic not clearly evidenced in code", 0.6))
  }
  if (has_filter_code) {
    return(partial("Filtering thresholds present in code but no saved object", 0.5))
  }
  fail("No evidence of cell filtering found")
}

grade_simple_04 <- function(work_dir) {
  content <- r_content(work_dir)
  if (grepl("NormalizeData", content) && grepl("10000|scale\\.factor", content)) {
    return(pass("Log-normalisation with scale factor found"))
  }
  if (grepl("NormalizeData", content)) {
    return(partial("NormalizeData used but scale factor unclear", 0.8))
  }
  fail("NormalizeData not found")
}

grade_simple_05 <- function(work_dir) {
  content <- r_content(work_dir)
  has_hvg <- grepl("FindVariableFeatures", content)
  has_plot <- grepl("VariableFeaturePlot|LabelPoints", content)
  if (has_hvg && has_plot) return(pass("HVG identified and visualised"))
  if (has_hvg) return(partial("HVG identified but no visualisation", 0.8))
  fail("FindVariableFeatures not found")
}

grade_simple_06 <- function(work_dir) {
  content <- r_content(work_dir)
  if (grepl("NormalizeData", content) && grepl("ScaleData", content)) {
    return(pass("NormalizeData + ScaleData found"))
  }
  if (grepl("ScaleData", content)) {
    return(partial("ScaleData found but no normalisation", 0.6))
  }
  fail("ScaleData not found")
}

grade_simple_07 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)
  has_pca_obj <- !is.null(seu) && "pca" %in% names(seu@reductions)
  has_elbow_code <- grepl("ElbowPlot", content)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if ((has_pca_obj || grepl("RunPCA", content)) && has_elbow_code && has_pdf) return(pass("PCA computed and visualised"))
  if (grepl("RunPCA", content) && has_elbow_code) return(partial("PCA + elbow in code but no output file confirmed", 0.8))
  if (has_pca_obj || grepl("RunPCA", content)) return(partial("PCA run but no elbow visualisation", 0.6))
  fail("RunPCA not found")
}

grade_simple_08 <- function(work_dir) {
  # No exact dataset-specific cluster-count range (was [5,15], a PBMC-3k
  # heuristic). Now checks for evidence of actual, non-trivial clustering:
  # a cluster assignment column that isn't a single uniform value.
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)

  if (!is.null(seu) && "seurat_clusters" %in% colnames(seu@meta.data)) {
    n_clusters <- length(unique(seu$seurat_clusters))
    if (n_clusters > 1) return(pass(sprintf("Clustering complete: %d clusters identified", n_clusters)))
    return(partial("Cluster column present but only a single cluster found", 0.5))
  }
  if (grepl("FindClusters", content)) {
    return(partial("FindClusters found in code but no saved/valid clustered object", 0.7))
  }
  fail("FindClusters not found")
}

grade_simple_09 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)
  has_umap_obj <- !is.null(seu) && "umap" %in% names(seu@reductions)
  has_umap_code <- grepl("RunUMAP", content) && grepl("DimPlot", content)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if ((has_umap_code && has_pdf) || (has_umap_obj && (grepl("DimPlot", content) || has_pdf))) {
    return(pass("UMAP computed and visualised"))
  }
  if (has_umap_code) return(partial("UMAP + DimPlot in code but no output file/object", 0.8))
  fail("RunUMAP not found")
}

grade_simple_10 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pdf <- length(find_pdf(work_dir)) > 0
  if (grepl("FeatureScatter", content) && has_pdf) return(pass("FeatureScatter plots generated"))
  if (grepl("FeatureScatter", content)) return(partial("FeatureScatter in code but no file", 0.8))
  fail("FeatureScatter not found")
}

grade_medium_01 <- function(work_dir) {
  # VERIFIED against prompts file: the prompt explicitly states "perform QC
  # filtering (200 < nFeature_RNA < 2500, percent.mt < 5)" — these are
  # task-specified thresholds, not a dataset assumption. Correctly kept exact.
  content <- r_content(work_dir)
  checks <- c(
    filter = grepl("200.*2500|nFeature_RNA.*2500", content) && grepl("percent\\.mt.*5", content),
    norm   = grepl("NormalizeData", content),
    hvg    = grepl("FindVariableFeatures", content),
    scale  = grepl("ScaleData", content)
  )
  n <- sum(checks)
  if (n == 4) return(pass("All 4 preprocessing steps present"))
  if (n >= 3) return(partial(sprintf("%d of 4 preprocessing steps present", n), 0.75))
  if (n >= 2) return(partial(sprintf("Only %d of 4 steps", n), 0.5))
  fail("Preprocessing steps not found")
}

grade_medium_02 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)
  has_pca <- (!is.null(seu) && "pca" %in% names(seu@reductions)) || grepl("RunPCA", content)
  has_umap <- (!is.null(seu) && "umap" %in% names(seu@reductions)) || grepl("RunUMAP", content)
  has_clusters <- (!is.null(seu) && "seurat_clusters" %in% colnames(seu@meta.data) && length(unique(seu$seurat_clusters)) > 1) ||
                  grepl("FindClusters", content)

  checks <- c(
    pca      = has_pca,
    knn      = grepl("FindNeighbors", content),
    cluster  = has_clusters,
    umap     = has_umap,
    dimplot  = grepl("DimPlot", content)
  )
  n <- sum(checks)
  if (n == 5) return(pass("Full dim reduction + clustering + UMAP pipeline"))
  if (n >= 3) return(partial(sprintf("%d of 5 steps correct", n), n / 5))
  fail("Core pipeline steps not found")
}

grade_medium_03 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_markers  <- grepl("FindAllMarkers", content)

  # No exact expected row count — just requires the CSV to actually have a
  # marker table structure with at least one real row.
  csv <- find_csv(work_dir)
  has_csv <- length(csv) > 0 && any(sapply(csv, function(f) {
    tryCatch({
      df <- read.csv(f)
      required_any <- c("gene", "Gene", "gene_symbol")
      has_gene_col <- any(required_any %in% colnames(df)) || any(grepl("gene", colnames(df), ignore.case = TRUE))
      has_cluster_col <- any(grepl("cluster", colnames(df), ignore.case = TRUE))
      has_effect_col <- any(grepl("log2?FC|logfc|avg_log", colnames(df), ignore.case = TRUE))
      has_gene_col && has_cluster_col && has_effect_col && nrow(df) > 0
    }, error = function(e) FALSE)
  }))

  if (has_pipeline && has_markers && has_csv) return(pass("Full pipeline + FindAllMarkers + valid marker CSV saved"))
  if (has_pipeline && has_markers) return(partial("Pipeline + FindAllMarkers found but no valid CSV", 0.8))
  if (has_pipeline) return(partial("Pipeline found but FindAllMarkers missing", 0.5))
  fail("Pipeline not found")
}

grade_medium_04 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_df <- grepl("doubletFinder|DoubletFinder", content)
  has_dynamic <- grepl("grep.*DF\\.classifications|DF\\.classifications.*grep", content)

  if (has_pipeline && has_df && has_dynamic) return(pass("Full pipeline + DoubletFinder with dynamic column detection"))
  if (has_pipeline && has_df) return(partial("Pipeline + DoubletFinder found (column may be hardcoded)", 0.8))
  if (has_df) return(partial("DoubletFinder found but pipeline incomplete", 0.5))
  fail("DoubletFinder not found")
}

grade_medium_05 <- function(work_dir) {
  # VERIFIED against prompts file: the prompt explicitly names "CD4 T cells
  # (IL7R, CCR7), CD14 Monocytes (CD14, LYZ), B cells (MS4A1), CD8 T cells
  # (CD8A), and NK cells (GNLY, NKG7)" — exact match to the marker list
  # below. Task-specified, correctly kept exact.
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_markers <- grepl("IL7R|CCR7|CD14|LYZ|MS4A1|CD8A|GNLY|NKG7", content)
  has_annotation <- grepl("RenameIdents|new\\.cluster\\.ids|cell.*type", content, ignore.case = TRUE)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if (has_pipeline && has_markers && has_annotation && has_pdf) return(pass("Full pipeline + annotation with canonical markers + UMAP"))
  if (has_pipeline && has_annotation) return(partial("Pipeline + annotation found", 0.8))
  if (has_pipeline && has_markers) return(partial("Pipeline + markers mentioned but no annotation", 0.6))
  fail("Cell type annotation not found")
}

grade_medium_06 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("FindClusters", content) && grepl("RunUMAP", content)
  has_markers <- grepl("FindAllMarkers", content)
  has_heatmap <- grepl("DoHeatmap", content)
  has_top10 <- grepl("slice_head.*10|head.*10|top.*10|n\\s*=\\s*10", content)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if (has_pipeline && has_markers && has_heatmap && has_top10 && has_pdf) return(pass("Full pipeline + top10 markers heatmap"))
  if (has_pipeline && has_heatmap) return(partial("Pipeline + DoHeatmap found", 0.8))
  fail("DoHeatmap not found")
}

grade_medium_07 <- function(work_dir) {
  # Cluster 0 vs cluster 3 is TASK-SPECIFIED (the prompt asks for this exact
  # comparison) — not a dataset-specific assumption, kept as-is.
  content <- r_content(work_dir)
  has_pipeline <- grepl("FindClusters", content)
  has_de <- grepl("FindMarkers", content) &&
            grepl("ident\\.1.*0|ident_1.*0", content) &&
            grepl("ident\\.2.*3|ident_2.*3", content)
  has_volcano <- grepl("volcano|ggplot|geom_point", content, ignore.case = TRUE)
  has_thresholds <- grepl("log2FC.*1|avg_log2FC.*1|p_adj|padj|p\\.val", content)

  if (has_pipeline && has_de && has_volcano && has_thresholds) return(pass("Cluster 0 vs 3 DE + volcano plot with thresholds"))
  if (has_pipeline && has_de) return(partial("Cluster 0 vs 3 DE found but no volcano", 0.7))
  if (has_de) return(partial("FindMarkers found but pipeline incomplete", 0.5))
  fail("Cluster 0 vs 3 FindMarkers not found")
}

grade_medium_08 <- function(work_dir) {
  # "pbmc_final.rds" naming is TASK-SPECIFIED (the prompt asks for this
  # exact filename). Reduction check kept as UMAP specifically — VERIFIED
  # against prompts file: the prompt explicitly says "clusters, computes
  # UMAP, and saves..." — UMAP is task-specified, not a dataset assumption.
  rds <- find_rds(work_dir)
  has_final_name <- length(rds) > 0 && any(grepl("pbmc_final|final", rds, ignore.case = TRUE))

  if (has_final_name) {
    result <- tryCatch({
      seu <- readRDS(rds[grepl("pbmc_final|final", rds, ignore.case = TRUE)][1])
      if (inherits(seu, "Seurat") && "umap" %in% names(seu@reductions)) {
        return(pass("pbmc_final.rds saved with a valid Seurat object and UMAP reduction"))
      }
      if (inherits(seu, "Seurat")) return(partial("pbmc_final.rds saved as valid Seurat object but no UMAP reduction found", 0.8))
      partial("pbmc_final.rds saved but could not verify contents", 0.7)
    }, error = function(e) partial("pbmc_final.rds saved but could not be loaded", 0.7))
    return(result)
  }
  if (length(rds) > 0) return(partial("RDS saved but not named pbmc_final.rds", 0.6))
  content <- r_content(work_dir)
  if (grepl("saveRDS|pbmc_final", content)) return(partial("saveRDS call found but file not written", 0.5))
  fail("No RDS saved")
}

grade_medium_09 <- function(work_dir) {
  # VERIFIED against prompts file: the prompt explicitly names "MS4A1,
  # GNLY, CD3E, CD14, FCER1A, FCGR3A, LYZ, PPBP, and CD8A" — exact match to
  # the marker list below. Task-specified, correctly kept exact.
  content <- r_content(work_dir)
  markers <- c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
  found <- sum(sapply(markers, function(m) grepl(m, content, fixed = TRUE)))
  has_pipeline <- grepl("FindClusters", content)
  has_feature <- grepl("FeaturePlot", content)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if (has_pipeline && has_feature && found >= 7 && has_pdf) return(pass(sprintf("FeaturePlot with %d/9 markers across full pipeline", found)))
  if (has_feature && found >= 5) return(partial(sprintf("FeaturePlot with %d/9 markers but pipeline incomplete", found), 0.7))
  fail("FeaturePlot or markers not found")
}

grade_medium_10 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("FindClusters", content)
  has_cc <- grepl("CellCycleScoring", content)
  has_cc_genes <- grepl("cc\\.genes|s\\.genes|g2m\\.genes", content)
  has_umap_cc <- grepl("Phase|group\\.by.*Phase", content)

  if (has_pipeline && has_cc && has_cc_genes && has_umap_cc) return(pass("Full pipeline + CellCycleScoring + UMAP by phase"))
  if (has_cc && has_cc_genes) return(partial("CellCycleScoring with gene sets found", 0.8))
  if (has_cc) return(partial("CellCycleScoring found but gene sets unclear", 0.6))
  fail("CellCycleScoring not found")
}

grade_hard_01 <- function(work_dir) {
  content <- r_content(work_dir)
  has_read10x <- grepl("Read10X", content)
  has_seurat <- grepl("CreateSeuratObject", content)
  has_pipeline <- grepl("FindClusters|RunPCA", content)

  if (has_read10x && has_seurat && has_pipeline) return(pass("Correctly inferred 10X format and ran pipeline"))
  if (has_read10x && has_seurat) return(partial("Correctly loaded data but pipeline incomplete", 0.7))
  if (has_read10x) return(partial("Read10X used but no Seurat object", 0.5))
  fail("Read10X not used — agent did not infer 10X format")
}

grade_hard_02 <- function(work_dir) {
  # No exact dataset-specific "cells remaining" bound (was <2700). Now
  # checks that QC + filtering happened and that the resulting object is
  # any valid, non-trivial size.
  seu <- load_seurat(work_dir)
  content <- r_content(work_dir)
  has_qc <- grepl("PercentageFeatureSet|percent\\.mt", content)
  has_filter <- grepl("subset.*nFeature|nFeature.*subset", content)

  if (has_qc && has_filter && !is.null(seu) && ncol(seu) > 0) {
    return(pass(sprintf("QC + filtering applied, %d cells remaining", ncol(seu))))
  }
  if (has_qc && has_filter) return(partial("QC + filtering code present but object not verifiable", 0.8))
  if (has_qc) return(partial("QC computed but no filtering", 0.5))
  fail("QC not performed")
}

grade_hard_03 <- function(work_dir) {
  content <- r_content(work_dir)
  checks <- c(
    norm  = grepl("NormalizeData|SCTransform", content),
    hvg   = grepl("FindVariableFeatures", content),
    scale = grepl("ScaleData", content),
    pca   = grepl("RunPCA", content),
    umap  = grepl("RunUMAP|RunTSNE", content)
  )
  n <- sum(checks)
  if (n == 5) return(pass("All dim reduction steps present"))
  if (n >= 3) return(partial(sprintf("%d of 5 dim reduction steps", n), n / 5))
  fail("Dim reduction steps not found")
}

grade_hard_04 <- function(work_dir) {
  # VERIFIED against prompts file: the prompt is genuinely vague — "I have
  # PBMC single-cell data and I want to identify different immune cell
  # populations. Write R code for the full analysis pipeline." — no cell
  # types are named. The check below is already correctly calibrated for
  # this: has_annotation is an OR across RenameIdents|SingleR|cell.*type|
  # CD4|CD8|B cell|Monocyte|NK — the named cell types are one of several
  # alternative signals, not a required match, so an agent using different
  # terminology still passes via the generic RenameIdents/cell.*type
  # branches. No change needed.
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_annotation <- grepl("RenameIdents|SingleR|cell.*type|CD4|CD8|B cell|Monocyte|NK", content, ignore.case = TRUE)
  has_umap <- grepl("RunUMAP|DimPlot", content)

  if (has_pipeline && has_annotation && has_umap) return(pass("Full pipeline + cell type identification"))
  if (has_pipeline && has_annotation) return(partial("Pipeline + annotation but no UMAP", 0.8))
  if (has_pipeline) return(partial("Pipeline found but no cell type annotation", 0.5))
  fail("Pipeline not found")
}

grade_hard_05 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pdf <- length(find_pdf(work_dir)) > 0
  if (grepl("FindClusters", content) && grepl("RunUMAP", content) && grepl("DimPlot", content) && has_pdf) {
    return(pass("Clustering + UMAP + visualisation"))
  }
  if (grepl("FindClusters", content) && grepl("RunUMAP", content)) return(partial("Cluster + UMAP in code", 0.8))
  fail("FindClusters or RunUMAP not found")
}

grade_hard_06 <- function(work_dir) {
  content <- r_content(work_dir)
  has_markers <- grepl("FindMarkers|FindAllMarkers", content)
  has_annotation <- grepl("RenameIdents|cell.*type|SingleR", content, ignore.case = TRUE)
  has_viz <- grepl("FeaturePlot|VlnPlot|DimPlot", content)

  if (has_markers && has_annotation && has_viz) return(pass("FindMarkers + annotation + visualisation"))
  if (has_markers && has_annotation) return(partial("Markers + annotation but no visualisation", 0.8))
  if (has_markers) return(partial("FindMarkers found but no annotation", 0.5))
  fail("FindMarkers not found")
}

grade_hard_07 <- function(work_dir) {
  content <- r_content(work_dir)
  checks <- c(
    load    = grepl("Read10X", content),
    qc      = grepl("PercentageFeatureSet", content),
    norm    = grepl("NormalizeData|SCTransform", content),
    pca     = grepl("RunPCA", content),
    cluster = grepl("FindClusters", content),
    umap    = grepl("RunUMAP", content),
    markers = grepl("FindAllMarkers|FindMarkers", content)
  )
  n <- sum(checks)
  if (n >= 6) return(pass(sprintf("Complete Seurat pipeline: %d of 7 steps", n)))
  if (n >= 4) return(partial(sprintf("%d of 7 pipeline steps", n), n / 7))
  fail(sprintf("Only %d of 7 pipeline steps found", n))
}

grade_hard_08 <- function(work_dir) {
  content <- r_content(work_dir)
  has_df <- grepl("doubletFinder|DoubletFinder", content)
  has_pipeline <- grepl("FindClusters", content)
  has_dynamic <- grepl("grep.*DF\\.classifications", content)
  has_subset <- grepl("Singlet|subset.*singlet", content, ignore.case = TRUE)

  if (has_pipeline && has_df && has_dynamic && has_subset) return(pass("Full pipeline + DoubletFinder + dynamic column + singlet subset"))
  if (has_df && has_subset) return(partial("DoubletFinder + subsetting found", 0.8))
  if (has_df) return(partial("DoubletFinder found but no subsetting", 0.5))
  fail("DoubletFinder not found")
}

grade_hard_09 <- function(work_dir) {
  content <- r_content(work_dir)
  has_de <- grepl("FindMarkers", content)
  has_viz <- grepl("ggplot|VlnPlot|FeaturePlot|volcano", content, ignore.case = TRUE)
  has_pdf <- length(find_pdf(work_dir)) > 0

  if (has_de && has_viz && has_pdf) return(pass("FindMarkers + visualisation + output file"))
  if (has_de && has_viz) return(partial("FindMarkers + visualisation in code", 0.8))
  if (has_de) return(partial("FindMarkers found but no visualisation", 0.5))
  fail("FindMarkers not found")
}

grade_hard_10 <- function(work_dir) {
  content <- r_content(work_dir)
  files <- list.files(work_dir, recursive = TRUE)
  exts <- tolower(tools::file_ext(files))

  checks <- c(
    pipeline = grepl("FindClusters", content) && grepl("RunUMAP", content),
    seed     = grepl("set\\.seed", content),
    save_rds = "rds" %in% exts,
    save_csv = "csv" %in% exts,
    save_pdf = any(c("pdf", "png") %in% exts),
    markers  = grepl("FindAllMarkers|FindMarkers", content)
  )
  n <- sum(checks)
  if (n >= 5) return(pass(sprintf("Production script: %d of 6 quality checks passed", n)))
  if (n >= 3) return(partial(sprintf("%d of 6 quality checks", n), n / 6))
  fail("Production script requirements not met")
}

# ── Dispatch table ───────────────────────────────────────────────────────────

GRADERS <- list(
  simple_01 = grade_simple_01, simple_02 = grade_simple_02, simple_03 = grade_simple_03,
  simple_04 = grade_simple_04, simple_05 = grade_simple_05, simple_06 = grade_simple_06,
  simple_07 = grade_simple_07, simple_08 = grade_simple_08, simple_09 = grade_simple_09,
  simple_10 = grade_simple_10,
  medium_01 = grade_medium_01, medium_02 = grade_medium_02, medium_03 = grade_medium_03,
  medium_04 = grade_medium_04, medium_05 = grade_medium_05, medium_06 = grade_medium_06,
  medium_07 = grade_medium_07, medium_08 = grade_medium_08, medium_09 = grade_medium_09,
  medium_10 = grade_medium_10,
  hard_01 = grade_hard_01, hard_02 = grade_hard_02, hard_03 = grade_hard_03,
  hard_04 = grade_hard_04, hard_05 = grade_hard_05, hard_06 = grade_hard_06,
  hard_07 = grade_hard_07, hard_08 = grade_hard_08, hard_09 = grade_hard_09,
  hard_10 = grade_hard_10
)

run_grader <- function(prompt_id, work_dir) {
  if (!dir.exists(work_dir)) {
    return(list(prompt_id = prompt_id, pass = FALSE, score = 0.0,
                reason = "work_dir does not exist",
                grounding_rate = 0.0, grounding_fns = list()))
  }

  grader_fn <- GRADERS[[prompt_id]]
  result <- if (is.null(grader_fn)) {
    list(pass = FALSE, score = 0.0, reason = paste("No grader defined for prompt_id:", prompt_id))
  } else {
    tryCatch(grader_fn(work_dir), error = function(e) {
      all_files <- tryCatch(list.files(work_dir, recursive = TRUE, full.names = FALSE), error = function(e2) character(0))
      list(pass = FALSE, score = 0.0,
           reason = paste("Grader error:", conditionMessage(e), "| Files in workdir:", paste(all_files, collapse = ", ")))
    })
  }

  grounding <- check_grounding(work_dir)

  list(
    prompt_id      = prompt_id,
    pass           = isTRUE(result$pass),
    score          = if (is.null(result$score)) 0.0 else as.numeric(result$score),
    reason         = if (is.null(result$reason) || result$reason == "") "No reason returned" else result$reason,
    grounding_rate = grounding$rate,
    grounding_fns  = grounding$found
  )
}

# ── CLI entry ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  prompt_id <- args[1]
  work_dir  <- normalizePath(args[2], winslash = "/", mustWork = FALSE)

  if (!dir.exists(work_dir)) {
    cat(toJSON(list(prompt_id = prompt_id, pass = FALSE, score = 0.0,
                     reason = paste("work_dir does not exist:", work_dir),
                     grounding_rate = 0.0, grounding_fns = list()), auto_unbox = TRUE), "\n")
    quit(status = 1)
  }

  result <- run_grader(prompt_id, work_dir)
  cat(toJSON(result, auto_unbox = TRUE), "\n")
}

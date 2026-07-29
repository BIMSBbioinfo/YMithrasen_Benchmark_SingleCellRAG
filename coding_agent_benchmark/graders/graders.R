# graders/graders.R
# Deterministic graders for coding agent benchmark
# Prompt IDs match internal benchmark: simple_01..10, medium_01..10, hard_01..10
# Each grader checks verifiable quantitative artifacts (scBench pattern)
# Called as: Rscript graders.R <prompt_id> <work_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})

# ── Shared helpers ─────────────────────────────────────────────────────────────

find_rds <- function(work_dir) {
  list.files(work_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
}

find_r_files <- function(work_dir) {
  files <- list.files(
    work_dir,
    pattern    = "\\.R$",
    full.names = TRUE,
    recursive  = TRUE
  )
  normalized <- gsub("\\\\", "/", files)
  files[!grepl("/recipes/", normalized, fixed = TRUE)]
}

find_csv <- function(work_dir) {
  list.files(work_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
}

find_pdf <- function(work_dir) {
  list.files(work_dir, pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)
}

r_content <- function(work_dir) {
  files <- find_r_files(work_dir)
  if (length(files) == 0) {
    # Debug: list ALL files in work_dir to help diagnose
    all_files <- list.files(work_dir, full.names = TRUE, recursive = TRUE, all.files = TRUE)
    return("")
  }
  paste(sapply(files, function(f) paste(readLines(f, warn=FALSE), collapse=" ")), collapse=" ")
}

load_seurat <- function(work_dir) {
  rds <- find_rds(work_dir)
  if (length(rds) == 0) return(NULL)
  for (f in rds) {
    tryCatch({
      obj <- readRDS(f)
      if (inherits(obj, "Seurat")) return(obj)
    }, error = function(e) NULL)
  }
  return(NULL)
}

pass  <- function(reason, score=1.0) list(pass=TRUE,  score=score, reason=reason)
fail  <- function(reason)            list(pass=FALSE, score=0.0,   reason=reason)
partial <- function(reason, score)   list(pass=FALSE, score=score, reason=reason)

# ── Grounding rate ─────────────────────────────────────────────────────────────

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
    if (nchar(content) == 0) return(list(rate=0.0, found=character(0)))
    found <- RECIPE_FUNCTIONS[sapply(RECIPE_FUNCTIONS, function(fn)
      grepl(fn, content, fixed=TRUE))]
    list(rate=length(found)/length(RECIPE_FUNCTIONS), found=found)
  }, error = function(e) {
    list(rate=0.0, found=character(0))
  })
}

# ── Individual graders ─────────────────────────────────────────────────────────

# simple_01: Load data + create Seurat object named pbmc
grade_simple_01 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  if (!is.null(seu) && ncol(seu) >= 2600 && nrow(seu) >= 13000)
    return(pass(paste("Seurat object: genes=", nrow(seu), "cells=", ncol(seu))))
  content <- r_content(work_dir)
  if (grepl("Read10X", content) && grepl("CreateSeuratObject", content))
    return(partial("Read10X + CreateSeuratObject used but no valid RDS found", 0.7))
  fail("Read10X and CreateSeuratObject not found")
}

# simple_02: QC metrics + violin plot
grade_simple_02 <- function(work_dir) {
  content <- r_content(work_dir)
  has_mt  <- grepl("PercentageFeatureSet|percent\\.mt", content)
  has_vln <- grepl("VlnPlot", content)
  has_pdf <- length(find_pdf(work_dir)) > 0
  if (has_mt && has_vln && has_pdf) return(pass("QC metrics and violin plot generated"))
  if (has_mt && has_vln)            return(partial("QC code correct but no plot file", 0.8))
  if (has_mt)                       return(partial("percent.mt computed only", 0.5))
  fail("QC metrics not found")
}

# simple_03: Load + filter (200 < nFeature < 2500, percent.mt < 5)
grade_simple_03 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  if (!is.null(seu)) {
    n <- ncol(seu)
    if (n >= 2000 && n <= 2700) return(pass(paste("Filtered cells:", n)))
    if (n < 2700)               return(partial(paste("Filtered but unusual count:", n), 0.7))
  }
  content <- r_content(work_dir)
  if (grepl("nFeature_RNA.*200|200.*nFeature_RNA", content) &&
      grepl("2500", content) && grepl("percent\\.mt.*5|5.*percent\\.mt", content))
    return(partial("Filter thresholds correct in code but no RDS", 0.7))
  fail("Filtering thresholds not found in code")
}

# simple_04: Load + log-normalise
grade_simple_04 <- function(work_dir) {
  content <- r_content(work_dir)
  if (grepl("NormalizeData", content) && grepl("10000|scale\\.factor", content))
    return(pass("Log-normalisation with scale.factor=10000 found"))
  if (grepl("NormalizeData", content))
    return(partial("NormalizeData used but scale factor unclear", 0.8))
  fail("NormalizeData not found")
}

# simple_05: HVG identification + visualisation
grade_simple_05 <- function(work_dir) {
  content <- r_content(work_dir)
  has_hvg  <- grepl("FindVariableFeatures", content) && grepl("2000|nfeatures", content)
  has_plot <- grepl("VariableFeaturePlot|LabelPoints", content)
  if (has_hvg && has_plot) return(pass("HVG identified and visualised"))
  if (has_hvg)             return(partial("HVG identified but no visualisation", 0.8))
  fail("FindVariableFeatures not found")
}

# simple_06: Load + normalise + scale
grade_simple_06 <- function(work_dir) {
  content <- r_content(work_dir)
  has_norm  <- grepl("NormalizeData", content)
  has_scale <- grepl("ScaleData",     content)
  if (has_norm && has_scale) return(pass("NormalizeData + ScaleData found"))
  if (has_scale)             return(partial("ScaleData found but no normalisation", 0.6))
  fail("ScaleData not found")
}

# simple_07: Full preprocessing + PCA + elbow plot
grade_simple_07 <- function(work_dir) {
  content  <- r_content(work_dir)
  has_pca  <- grepl("RunPCA",    content)
  has_elbow<- grepl("ElbowPlot", content)
  has_pdf  <- length(find_pdf(work_dir)) > 0
  if (has_pca && has_elbow && has_pdf) return(pass("PCA + elbow plot generated"))
  if (has_pca && has_elbow)            return(partial("PCA + elbow in code but no file", 0.8))
  if (has_pca)                         return(partial("PCA run but no elbow plot", 0.6))
  fail("RunPCA not found")
}

# simple_08: Full preprocessing + clustering (resolution=0.5, 10 PCs)
grade_simple_08 <- function(work_dir) {
  seu <- load_seurat(work_dir)
  if (!is.null(seu) && "seurat_clusters" %in% colnames(seu@meta.data)) {
    n_clusters <- length(levels(seu$seurat_clusters))
    if (n_clusters >= 5 && n_clusters <= 15)
      return(pass(paste("Clustering complete:", n_clusters, "clusters")))
  }
  content <- r_content(work_dir)
  if (grepl("FindClusters", content) && grepl("0\\.5", content))
    return(partial("FindClusters(resolution=0.5) found but no RDS", 0.7))
  fail("FindClusters not found")
}

# simple_09: Full pipeline + UMAP + cluster visualisation
grade_simple_09 <- function(work_dir) {
  content  <- r_content(work_dir)
  has_umap <- grepl("RunUMAP",  content)
  has_plot <- grepl("DimPlot",  content)
  has_pdf  <- length(find_pdf(work_dir)) > 0
  seu <- load_seurat(work_dir)
  has_umap_in_obj <- !is.null(seu) && "umap" %in% names(seu@reductions)
  if ((has_umap && has_plot && has_pdf) || has_umap_in_obj)
    return(pass("UMAP computed and visualised"))
  if (has_umap && has_plot)
    return(partial("UMAP + DimPlot in code but no output file", 0.8))
  fail("RunUMAP not found")
}

# simple_10: QC scatter plots
grade_simple_10 <- function(work_dir) {
  content  <- r_content(work_dir)
  has_scatter <- grepl("FeatureScatter", content)
  has_pdf     <- length(find_pdf(work_dir)) > 0
  if (has_scatter && has_pdf) return(pass("FeatureScatter plots generated"))
  if (has_scatter)            return(partial("FeatureScatter in code but no file", 0.8))
  fail("FeatureScatter not found")
}

# medium_01: QC filter + normalise + HVG + scale
grade_medium_01 <- function(work_dir) {
  content <- r_content(work_dir)
  checks  <- c(
    filter  = grepl("200.*2500|nFeature_RNA.*2500", content) && grepl("percent\\.mt.*5", content),
    norm    = grepl("NormalizeData", content),
    hvg     = grepl("FindVariableFeatures", content),
    scale   = grepl("ScaleData", content)
  )
  n <- sum(checks)
  if (n == 4) return(pass("All 4 preprocessing steps present"))
  if (n >= 3) return(partial(paste("Steps present:", paste(names(checks[checks]), collapse=",")), 0.75))
  if (n >= 2) return(partial(paste("Only", n, "of 4 steps"), 0.5))
  fail("Preprocessing steps not found")
}

# medium_02: preprocessing + PCA + KNN + cluster (10 PCs, res 0.5) + UMAP
grade_medium_02 <- function(work_dir) {
  content <- r_content(work_dir)
  checks <- c(
    pca      = grepl("RunPCA",       content),
    knn      = grepl("FindNeighbors",content),
    cluster  = grepl("FindClusters", content) && grepl("0\\.5", content),
    umap     = grepl("RunUMAP",      content),
    dimplot  = grepl("DimPlot",      content)
  )
  n <- sum(checks)
  if (n == 5) return(pass("Full dim reduction + clustering + UMAP pipeline"))
  if (n >= 3) return(partial(paste(n, "of 5 steps correct"), n/5))
  fail("Core pipeline steps not found")
}

# medium_03: Full pipeline + FindAllMarkers (logfc=0.25)
grade_medium_03 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_markers  <- grepl("FindAllMarkers", content) && grepl("0\\.25", content)
  csv <- find_csv(work_dir)
  has_csv <- length(csv) > 0 && any(sapply(csv, function(f) {
    tryCatch({
      df <- read.csv(f)
      all(c("gene","cluster","avg_log2FC") %in% colnames(df)) && nrow(df) > 50
    }, error=function(e) FALSE)
  }))
  if (has_pipeline && has_markers && has_csv)
    return(pass("Full pipeline + FindAllMarkers + CSV saved"))
  if (has_pipeline && has_markers)
    return(partial("Pipeline + FindAllMarkers found but no CSV", 0.8))
  if (has_pipeline)
    return(partial("Pipeline found but FindAllMarkers missing", 0.5))
  fail("Pipeline not found")
}

# medium_04: Full pipeline + DoubletFinder
grade_medium_04 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_df       <- grepl("doubletFinder|DoubletFinder", content)
  has_dynamic  <- grepl("grep.*DF\\.classifications|DF\\.classifications.*grep", content)
  if (has_pipeline && has_df && has_dynamic)
    return(pass("Full pipeline + DoubletFinder with dynamic column detection"))
  if (has_pipeline && has_df)
    return(partial("Pipeline + DoubletFinder found (column may be hardcoded)", 0.8))
  if (has_df)
    return(partial("DoubletFinder found but pipeline incomplete", 0.5))
  fail("DoubletFinder not found")
}

# medium_05: Full pipeline + manual cell type annotation + annotated UMAP
grade_medium_05 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline  <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_markers   <- grepl("IL7R|CCR7|CD14|LYZ|MS4A1|CD8A|GNLY|NKG7", content)
  has_annotation<- grepl("RenameIdents|new\\.cluster\\.ids|cell.*type", content, ignore.case=TRUE)
  has_pdf       <- length(find_pdf(work_dir)) > 0
  if (has_pipeline && has_markers && has_annotation && has_pdf)
    return(pass("Full pipeline + annotation with canonical markers + UMAP"))
  if (has_pipeline && has_annotation)
    return(partial("Pipeline + annotation found", 0.8))
  if (has_pipeline && has_markers)
    return(partial("Pipeline + markers mentioned but no annotation", 0.6))
  fail("Cell type annotation not found")
}

# medium_06: Full pipeline + FindAllMarkers + DoHeatmap top 10
grade_medium_06 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("FindClusters", content) && grepl("RunUMAP", content)
  has_markers  <- grepl("FindAllMarkers", content)
  has_heatmap  <- grepl("DoHeatmap", content)
  has_top10    <- grepl("slice_head.*10|head.*10|top.*10|n.*=.*10", content)
  has_pdf      <- length(find_pdf(work_dir)) > 0
  if (has_pipeline && has_markers && has_heatmap && has_top10 && has_pdf)
    return(pass("Full pipeline + top10 markers heatmap"))
  if (has_pipeline && has_heatmap)
    return(partial("Pipeline + DoHeatmap found", 0.8))
  fail("DoHeatmap not found")
}

# medium_07: Full pipeline + FindMarkers cluster0 vs cluster3 + volcano plot
grade_medium_07 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline <- grepl("FindClusters", content)
  has_de       <- grepl("FindMarkers", content) &&
                  grepl("ident\\.1.*0|ident_1.*0", content) &&
                  grepl("ident\\.2.*3|ident_2.*3", content)
  has_volcano  <- grepl("volcano|ggplot|geom_point", content, ignore.case=TRUE)
  has_thresholds <- grepl("log2FC.*1|avg_log2FC.*1|p_adj|padj|p\\.val", content)
  if (has_pipeline && has_de && has_volcano && has_thresholds)
    return(pass("Cluster 0 vs 3 DE + volcano plot with thresholds"))
  if (has_pipeline && has_de)
    return(partial("Cluster 0 vs 3 DE found but no volcano", 0.7))
  if (has_de) return(partial("FindMarkers found but pipeline incomplete", 0.5))
  fail("Cluster 0 vs 3 FindMarkers not found")
}

# medium_08: Complete pipeline + save pbmc_final.rds
grade_medium_08 <- function(work_dir) {
  rds <- find_rds(work_dir)
  has_final <- any(grepl("pbmc_final|final", rds, ignore.case=TRUE))
  if (has_final) {
    tryCatch({
      seu <- readRDS(rds[grepl("final", rds, ignore.case=TRUE)][1])
      if (inherits(seu, "Seurat") && "umap" %in% names(seu@reductions))
        return(pass(paste("pbmc_final.rds saved with UMAP, cells:", ncol(seu))))
    }, error=function(e) NULL)
    return(partial("pbmc_final.rds saved but could not verify contents", 0.8))
  }
  if (length(rds) > 0) return(partial("RDS saved but not named pbmc_final.rds", 0.7))
  content <- r_content(work_dir)
  if (grepl("saveRDS|pbmc_final", content))
    return(partial("saveRDS call found but file not written", 0.5))
  fail("No RDS saved")
}

# medium_09: Full pipeline + FeaturePlot for 9 canonical PBMC markers
grade_medium_09 <- function(work_dir) {
  content <- r_content(work_dir)
  markers <- c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A","LYZ","PPBP","CD8A")
  has_pipeline  <- grepl("FindClusters", content)
  has_feature   <- grepl("FeaturePlot",  content)
  markers_found <- sum(sapply(markers, function(m) grepl(m, content)))
  has_pdf       <- length(find_pdf(work_dir)) > 0
  if (has_pipeline && has_feature && markers_found >= 7 && has_pdf)
    return(pass(paste("FeaturePlot with", markers_found, "of 9 markers")))
  if (has_feature && markers_found >= 5)
    return(partial(paste("FeaturePlot with", markers_found, "markers but pipeline incomplete"), 0.7))
  fail(paste("FeaturePlot or markers not found (found", markers_found, "markers)"))
}

# medium_10: Full pipeline + CellCycleScoring + UMAP coloured by phase
grade_medium_10 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline  <- grepl("FindClusters", content)
  has_cc        <- grepl("CellCycleScoring", content)
  has_cc_genes  <- grepl("cc\\.genes|s\\.genes|g2m\\.genes", content)
  has_umap_cc   <- grepl("Phase|group\\.by.*Phase", content)
  if (has_pipeline && has_cc && has_cc_genes && has_umap_cc)
    return(pass("Full pipeline + CellCycleScoring + UMAP by phase"))
  if (has_cc && has_cc_genes)
    return(partial("CellCycleScoring with gene sets found", 0.8))
  if (has_cc)
    return(partial("CellCycleScoring found but gene sets unclear", 0.6))
  fail("CellCycleScoring not found")
}

# hard prompts: vague/underspecified — test whether agent resolves ambiguity correctly

# hard_01: "3 files in working directory" — must infer 10X format and use Read10X
grade_hard_01 <- function(work_dir) {
  content <- r_content(work_dir)
  has_read10x  <- grepl("Read10X",         content)
  has_seurat   <- grepl("CreateSeuratObject", content)
  has_pipeline <- grepl("FindClusters|RunPCA", content)
  if (has_read10x && has_seurat && has_pipeline)
    return(pass("Correctly inferred 10X format and ran pipeline"))
  if (has_read10x && has_seurat)
    return(partial("Correctly loaded data but pipeline incomplete", 0.7))
  if (has_read10x)
    return(partial("Read10X used but no Seurat object", 0.5))
  fail("Read10X not used — agent did not infer 10X format")
}

# hard_02: "QC and remove low-quality cells" — vague thresholds
grade_hard_02 <- function(work_dir) {
  content <- r_content(work_dir)
  has_qc     <- grepl("PercentageFeatureSet|percent\\.mt", content)
  has_filter <- grepl("subset.*nFeature|nFeature.*subset", content)
  seu <- load_seurat(work_dir)
  filtered_correctly <- !is.null(seu) && ncol(seu) < 2700
  if (has_qc && has_filter && filtered_correctly)
    return(pass(paste("QC + filtering applied, cells remaining:", if(!is.null(seu)) ncol(seu) else "?")))
  if (has_qc && has_filter)
    return(partial("QC + filtering code present", 0.8))
  if (has_qc)
    return(partial("QC computed but no filtering", 0.5))
  fail("QC not performed")
}

# hard_03: "normalise, HVG, dimensionality reduction" — must choose appropriate methods
grade_hard_03 <- function(work_dir) {
  content <- r_content(work_dir)
  checks <- c(
    norm   = grepl("NormalizeData|SCTransform", content),
    hvg    = grepl("FindVariableFeatures",      content),
    scale  = grepl("ScaleData",                 content),
    pca    = grepl("RunPCA",                    content),
    umap   = grepl("RunUMAP|RunTSNE",           content)
  )
  n <- sum(checks)
  if (n == 5) return(pass("All dim reduction steps present"))
  if (n >= 3) return(partial(paste(n, "of 5 dim reduction steps"), n/5))
  fail("Dim reduction steps not found")
}

# hard_04: "identify different immune cell populations" — full pipeline + annotation
grade_hard_04 <- function(work_dir) {
  content <- r_content(work_dir)
  has_pipeline   <- grepl("RunPCA", content) && grepl("FindClusters", content)
  has_annotation <- grepl("RenameIdents|SingleR|cell.*type|CD4|CD8|B cell|Monocyte|NK",
                           content, ignore.case=TRUE)
  has_umap       <- grepl("RunUMAP|DimPlot", content)
  if (has_pipeline && has_annotation && has_umap)
    return(pass("Full pipeline + cell type identification"))
  if (has_pipeline && has_annotation)
    return(partial("Pipeline + annotation but no UMAP", 0.8))
  if (has_pipeline)
    return(partial("Pipeline found but no cell type annotation", 0.5))
  fail("Pipeline not found")
}

# hard_05: "cluster and visualise on UMAP" — core workflow
grade_hard_05 <- function(work_dir) {
  content <- r_content(work_dir)
  has_cluster <- grepl("FindClusters", content)
  has_umap    <- grepl("RunUMAP",      content)
  has_plot    <- grepl("DimPlot",      content)
  has_pdf     <- length(find_pdf(work_dir)) > 0
  if (has_cluster && has_umap && has_plot && has_pdf)
    return(pass("Clustering + UMAP + visualisation"))
  if (has_cluster && has_umap)
    return(partial("Cluster + UMAP in code", 0.8))
  fail("FindClusters or RunUMAP not found")
}

# hard_06: "find marker genes and annotate cell types"
grade_hard_06 <- function(work_dir) {
  content <- r_content(work_dir)
  has_markers    <- grepl("FindMarkers|FindAllMarkers", content)
  has_annotation <- grepl("RenameIdents|cell.*type|SingleR", content, ignore.case=TRUE)
  has_viz        <- grepl("FeaturePlot|VlnPlot|DimPlot", content)
  if (has_markers && has_annotation && has_viz)
    return(pass("FindMarkers + annotation + visualisation"))
  if (has_markers && has_annotation)
    return(partial("Markers + annotation but no visualisation", 0.8))
  if (has_markers)
    return(partial("FindMarkers found but no annotation", 0.5))
  fail("FindMarkers not found")
}

# hard_07: "complete R analysis pipeline using Seurat"
grade_hard_07 <- function(work_dir) {
  content <- r_content(work_dir)
  steps <- c(
    load    = grepl("Read10X",             content),
    qc      = grepl("PercentageFeatureSet",content),
    norm    = grepl("NormalizeData|SCTransform", content),
    pca     = grepl("RunPCA",              content),
    cluster = grepl("FindClusters",        content),
    umap    = grepl("RunUMAP",             content),
    markers = grepl("FindAllMarkers|FindMarkers", content)
  )
  n <- sum(steps)
  if (n >= 6) return(pass(paste("Complete Seurat pipeline:", n, "of 7 steps")))
  if (n >= 4) return(partial(paste(n, "of 7 pipeline steps"), n/7))
  fail(paste("Only", n, "of 7 pipeline steps found"))
}

# hard_08: "detect and remove doublets"
grade_hard_08 <- function(work_dir) {
  content <- r_content(work_dir)
  has_df       <- grepl("doubletFinder|DoubletFinder", content)
  has_pipeline <- grepl("FindClusters", content)
  has_dynamic  <- grepl("grep.*DF\\.classifications", content)
  has_subset   <- grepl("Singlet|subset.*singlet", content, ignore.case=TRUE)
  if (has_pipeline && has_df && has_dynamic && has_subset)
    return(pass("Full pipeline + DoubletFinder + dynamic column + singlet subset"))
  if (has_df && has_subset)
    return(partial("DoubletFinder + subsetting found", 0.8))
  if (has_df)
    return(partial("DoubletFinder found but no subsetting", 0.5))
  fail("DoubletFinder not found")
}

# hard_09: "differential expression between two clusters + visualise"
grade_hard_09 <- function(work_dir) {
  content <- r_content(work_dir)
  has_de   <- grepl("FindMarkers", content)
  has_viz  <- grepl("ggplot|VlnPlot|FeaturePlot|volcano", content, ignore.case=TRUE)
  has_pdf  <- length(find_pdf(work_dir)) > 0
  if (has_de && has_viz && has_pdf) return(pass("FindMarkers + visualisation + output file"))
  if (has_de && has_viz)            return(partial("FindMarkers + visualisation in code", 0.8))
  if (has_de)                       return(partial("FindMarkers found but no visualisation", 0.5))
  fail("FindMarkers not found")
}

# hard_10: "production-ready script + save all outputs"
grade_hard_10 <- function(work_dir) {
  content <- r_content(work_dir)
  rds <- find_rds(work_dir)
  csv <- find_csv(work_dir)
  pdf <- find_pdf(work_dir)
  steps <- c(
    pipeline  = grepl("FindClusters", content) && grepl("RunUMAP", content),
    seed      = grepl("set\\.seed", content),
    save_rds  = length(rds) > 0,
    save_csv  = length(csv) > 0,
    save_pdf  = length(pdf) > 0,
    markers   = grepl("FindAllMarkers|FindMarkers", content)
  )
  n <- sum(steps)
  if (n >= 5) return(pass(paste("Production script:", n, "of 6 quality checks passed")))
  if (n >= 3) return(partial(paste(n, "of 6 quality checks"), n/6))
  fail("Production script requirements not met")
}

# ── Dispatcher ─────────────────────────────────────────────────────────────────

run_grader <- function(prompt_id, work_dir) {
  # Safety: ensure work_dir exists
  if (!dir.exists(work_dir)) {
    return(list(prompt_id=prompt_id, pass=FALSE, score=0.0,
                reason="work_dir does not exist",
                grounding_rate=0.0, grounding_fns=list()))
  }
  grade_fn <- switch(prompt_id,
    "simple_01" = grade_simple_01, "simple_02" = grade_simple_02,
    "simple_03" = grade_simple_03, "simple_04" = grade_simple_04,
    "simple_05" = grade_simple_05, "simple_06" = grade_simple_06,
    "simple_07" = grade_simple_07, "simple_08" = grade_simple_08,
    "simple_09" = grade_simple_09, "simple_10" = grade_simple_10,
    "medium_01" = grade_medium_01, "medium_02" = grade_medium_02,
    "medium_03" = grade_medium_03, "medium_04" = grade_medium_04,
    "medium_05" = grade_medium_05, "medium_06" = grade_medium_06,
    "medium_07" = grade_medium_07, "medium_08" = grade_medium_08,
    "medium_09" = grade_medium_09, "medium_10" = grade_medium_10,
    "hard_01"   = grade_hard_01,   "hard_02"   = grade_hard_02,
    "hard_03"   = grade_hard_03,   "hard_04"   = grade_hard_04,
    "hard_05"   = grade_hard_05,   "hard_06"   = grade_hard_06,
    "hard_07"   = grade_hard_07,   "hard_08"   = grade_hard_08,
    "hard_09"   = grade_hard_09,   "hard_10"   = grade_hard_10,
    stop(paste("Unknown prompt_id:", prompt_id))
  )
  result <- tryCatch(grade_fn(work_dir), error = function(e) {
    # On grader error, check if R files exist and report debug info
    all_files <- tryCatch(
      list.files(work_dir, recursive=TRUE, full.names=FALSE),
      error = function(e2) character(0)
    )
    list(pass=FALSE, score=0.0,
         reason=paste("Grader error:", conditionMessage(e),
                      "| Files in workdir:", paste(all_files, collapse=", ")))
  })
  grounding <- check_grounding(work_dir)
  list(
    prompt_id      = prompt_id,
    pass           = isTRUE(result$pass),
    score          = if (is.null(result$score)) 0.0 else as.numeric(result$score),
    reason         = if (is.null(result$reason) || result$reason == "") {
                       # Report what files exist when reason is empty
                       all_f <- tryCatch(
                         list.files(work_dir, recursive=TRUE, full.names=FALSE),
                         error=function(e) character(0))
                       paste("No reason returned. Files:", paste(all_f, collapse=", "))
                     } else result$reason,
    grounding_rate = grounding$rate,
    grounding_fns  = grounding$found
  )
}

# ── CLI entry ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  prompt_id <- args[1]
  # Normalize path separators for R on Windows
  work_dir  <- normalizePath(args[2], winslash = "/", mustWork = FALSE)

  if (!dir.exists(work_dir)) {
    cat(toJSON(list(
      prompt_id      = prompt_id,
      pass           = FALSE,
      score          = 0.0,
      reason         = paste("work_dir does not exist:", work_dir),
      grounding_rate = 0.0,
      grounding_fns  = list()
    ), auto_unbox = TRUE), "\n")
    quit(status = 1)
  }

  result <- run_grader(prompt_id, work_dir)
  cat(toJSON(result, auto_unbox = TRUE), "\n")
}
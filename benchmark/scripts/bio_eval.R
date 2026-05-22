# Biological evaluation metrics for the ROLV benchmark.
# Run AFTER the LLM has generated and executed code that produces a Seurat object.
#
# Usage:
#   Rscript bio_eval.R \
#     --rds_path /path/to/output_seurat.rds \
#     --gt_labels /path/to/ground_truth.csv \
#     --out_json  /path/to/bio_metrics.json \
#     --run_batch TRUE \
#     --batch_col orig.ident
#
# Outputs a JSON file with all biological metrics for import into benchmark_runner.js

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(jsonlite)
})

#  CLI args ─
option_list <- list(
  make_option("--rds_path",  type="character", help="Path to output Seurat .rds file"),
  make_option("--gt_labels", type="character", help="CSV with columns: cell_barcode, cell_type"),
  make_option("--out_json",  type="character", default="bio_metrics.json"),
  make_option("--run_batch", type="logical",   default=FALSE, help="Run batch correction metrics"),
  make_option("--batch_col", type="character", default="orig.ident", help="Batch column in Seurat metadata"),
  make_option("--cluster_col", type="character", default="seurat_clusters"),
  make_option("--embedding_col", type="character", default="umap")
)
parser <- OptionParser(option_list=option_list)
opts   <- parse_args(parser)

if (is.null(opts$rds_path)) {
  cat('Usage: Rscript bio_eval.R --rds_path <file.rds> --gt_labels <labels.csv>\n')
  quit(status=1)
}

results <- list()

#  Load Seurat object
cat("[bio_eval] Loading Seurat object:", opts$rds_path, "\n")
seurat_obj <- tryCatch(
  readRDS(opts$rds_path),
  error = function(e) { cat("[ERROR] Cannot load .rds:", e$message, "\n"); quit(status=1) }
)
cat("[bio_eval] Cells:", ncol(seurat_obj), " | Features:", nrow(seurat_obj), "\n")

#  METRIC 1: Clustering Quality — ARI + F1
if (!is.null(opts$gt_labels) && file.exists(opts$gt_labels)) {
  cat("\n[Metric 1] Clustering Quality (ARI, F1)\n")

  gt <- read.csv(opts$gt_labels, row.names=1, stringsAsFactors=FALSE)
  # Align barcodes
  shared_cells <- intersect(rownames(gt), colnames(seurat_obj))
  cat("  Shared cells:", length(shared_cells), "\n")

  if (length(shared_cells) > 10) {
    pred_clusters <- as.character(seurat_obj@meta.data[shared_cells, opts$cluster_col])
    true_labels   <- gt[shared_cells, "cell_type"]

    # Install mclust if needed (provides adjustedRandIndex)
    if (!requireNamespace("mclust", quietly=TRUE)) {
      install.packages("mclust", repos="https://cloud.r-project.org", quiet=TRUE)
    }
    suppressPackageStartupMessages(library(mclust))
    ari <- adjustedRandIndex(pred_clusters, true_labels)
    cat("  ARI:", round(ari, 4), "\n")

    # Hungarian-matched F1
    # Build contingency table and use linear_sum_assignment equivalent
    cont_table <- table(pred=pred_clusters, true=true_labels)
    # Greedy matching: assign each predicted cluster to best matching true label
    matched_f1_values <- c()
    assigned_true <- c()
    for (pred_cl in rownames(cont_table)) {
      row <- cont_table[pred_cl,]
      row[assigned_true] <- 0  # don't reuse
      best_true <- names(which.max(row))
      tp <- cont_table[pred_cl, best_true]
      fp <- sum(cont_table[pred_cl,]) - tp
      fn <- sum(cont_table[, best_true]) - tp
      prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
      rec  <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
      f1   <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
      matched_f1_values <- c(matched_f1_values, f1)
      assigned_true <- c(assigned_true, best_true)
    }
    macro_f1 <- mean(matched_f1_values)
    cat("  Macro F1 (matched):", round(macro_f1, 4), "\n")

    results$clustering <- list(
      ari       = round(ari, 4),
      macro_f1  = round(macro_f1, 4),
      n_cells   = length(shared_cells),
      n_clusters = length(unique(pred_clusters))
    )
  } else {
    cat("  [SKIP] Not enough shared cells\n")
  }
}

#  METRIC 2: Marker Gene Agreement
cat("\n[Metric 2] Marker Gene Agreement\n")

# CellMarker 2.0 — a small curated subset for PBMC cell types
# For full analysis, download CellMarker2 from: http://bio-bigdata.hrbmu.edu.cn/CellMarker/
PBMC_MARKERS <- list(
  "CD4 T cell"     = c("IL7R","CCR7","CD4","TCF7","LEF1"),
  "CD8 T cell"     = c("CD8A","CD8B","GZMB","GZMA","PRF1"),
  "B cell"         = c("MS4A1","CD79A","CD79B","BANK1","IGHM"),
  "NK cell"        = c("GNLY","NKG7","NCAM1","KLRD1","FGFBP2"),
  "CD14 Monocyte"  = c("CD14","LYZ","S100A8","S100A9","FCGR3A"),
  "Dendritic cell" = c("FCER1A","CST3","ITGAX","CLEC9A","XCR1"),
  "Platelet"       = c("PPBP","PF4","GP1BA","ITGA2B"),
  "T reg"          = c("FOXP3","IL2RA","CTLA4","TIGIT")
)

marker_scores <- list()
tryCatch({
  # Find all markers from the Seurat object
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[opts$cluster_col]]
  all_markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos    = TRUE,
    min.pct     = 0.25,
    logfc.threshold = 0.25,
    verbose     = FALSE
  )

  # For each cluster, get top 20 markers and compute overlap with each cell type
  for (cl in unique(all_markers$cluster)) {
    top_genes <- all_markers[all_markers$cluster == cl, ]
    top_genes <- head(top_genes[order(top_genes$avg_log2FC, decreasing=TRUE), "gene"], 20)
    best_type  <- ""
    best_score <- 0
    for (ct in names(PBMC_MARKERS)) {
      overlap <- length(intersect(top_genes, PBMC_MARKERS[[ct]]))
      score   <- overlap / length(PBMC_MARKERS[[ct]])
      if (score > best_score) { best_score <- score; best_type <- ct }
    }
    marker_scores[[as.character(cl)]] <- list(
      best_match   = best_type,
      overlap_score = round(best_score, 4),
      top_genes    = top_genes[1:min(5, length(top_genes))]
    )
    cat(sprintf("  Cluster %s → %s (overlap=%.2f)\n", cl, best_type, best_score))
  }

  avg_marker_overlap <- mean(sapply(marker_scores, function(x) x$overlap_score))
  cat("  Average marker overlap:", round(avg_marker_overlap, 4), "\n")
  results$marker_gene <- list(per_cluster=marker_scores, average_overlap=round(avg_marker_overlap, 4))

}, error=function(e) {
  cat("  [WARN] Marker gene evaluation failed:", e$message, "\n")
})

#  METRIC 3: Batch Correction — kBET, LISI, ASW
if (opts$run_batch) {
  cat("\n[Metric 3] Batch Correction Metrics\n")

  batch_meta <- seurat_obj@meta.data[[opts$batch_col]]
  if (is.null(batch_meta)) {
    cat("  [SKIP] batch_col", opts$batch_col, "not found in metadata\n")
  } else {

    # Get embedding
    embed_key <- toupper(opts$embedding_col)
    embedding <- tryCatch(
      Embeddings(seurat_obj, reduction=opts$embedding_col)[, 1:2],
      error=function(e) NULL
    )
    if (is.null(embedding)) {
      embedding <- tryCatch(
        Embeddings(seurat_obj, reduction="pca")[, 1:10],
        error=function(e) NULL
      )
      cat("  [INFO] Using PCA embedding for batch metrics\n")
    }

    batch_results <- list()

    #  3a. kBET
    tryCatch({
      if (!requireNamespace("kBET", quietly=TRUE)) {
        # Install from GitHub — requires remotes
        if (!requireNamespace("remotes", quietly=TRUE))
          install.packages("remotes", repos="https://cloud.r-project.org", quiet=TRUE)
        remotes::install_github("theislab/kBET", quiet=TRUE)
      }
      suppressPackageStartupMessages(library(kBET))
      batch_factor <- as.factor(batch_meta)
      kbet_result  <- kBET(embedding, batch=batch_factor, plot=FALSE, verbose=FALSE)
      kbet_score   <- 1 - kbet_result$summary$kBET.observed[1]  # acceptance rate
      cat("  kBET acceptance rate:", round(kbet_score, 4), "\n")
      batch_results$kBET <- round(kbet_score, 4)
    }, error=function(e) cat("  [WARN] kBET failed:", e$message, "\n"))

    #  3b. LISI (iLISI + cLISI) ─
    tryCatch({
      if (!requireNamespace("lisi", quietly=TRUE)) {
        if (!requireNamespace("remotes", quietly=TRUE))
          install.packages("remotes", repos="https://cloud.r-project.org", quiet=TRUE)
        remotes::install_github("immunogenomics/lisi", quiet=TRUE)
      }
      suppressPackageStartupMessages(library(lisi))

      meta_df <- data.frame(
        batch     = batch_meta,
        cell_type = seurat_obj@meta.data[[opts$cluster_col]]
      )
      lisi_scores <- compute_lisi(embedding, meta_df, c("batch", "cell_type"))
      iLISI <- mean(lisi_scores$batch)      # higher = better batch mixing
      cLISI <- mean(lisi_scores$cell_type)  # lower = better cell type separation
      cat("  iLISI:", round(iLISI, 4), "(higher=better batch mixing)\n")
      cat("  cLISI:", round(cLISI, 4), "(lower=better cell type separation)\n")
      batch_results$iLISI <- round(iLISI, 4)
      batch_results$cLISI <- round(cLISI, 4)
    }, error=function(e) cat("  [WARN] LISI failed:", e$message, "\n"))

    #  3c. ASW (Average Silhouette Width)
    tryCatch({
      suppressPackageStartupMessages(library(cluster))
      cluster_labels <- as.integer(as.factor(seurat_obj@meta.data[[opts$cluster_col]]))
      # Sample max 5000 cells for speed
      n <- nrow(embedding)
      idx <- if (n > 5000) sample(n, 5000) else seq_len(n)
      dist_mat <- dist(embedding[idx, ])
      sil      <- silhouette(cluster_labels[idx], dist_mat)
      asw      <- mean(sil[, 3])
      cat("  ASW:", round(asw, 4), "(higher=better separation)\n")
      batch_results$ASW <- round(asw, 4)
    }, error=function(e) cat("  [WARN] ASW failed:", e$message, "\n"))

    results$batch_correction <- batch_results
  }
}

#  WRITE OUTPUT JSON
results$n_cells    <- ncol(seurat_obj)
results$n_features <- nrow(seurat_obj)
results$evaluated_at <- Sys.time()

out <- toJSON(results, auto_unbox=TRUE, pretty=TRUE)
writeLines(out, opts$out_json)
cat("\n[bio_eval] Results written →", opts$out_json, "\n")

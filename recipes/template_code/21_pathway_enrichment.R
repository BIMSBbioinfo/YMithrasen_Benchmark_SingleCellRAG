library(Seurat)
library(enrichR)
library(dplyr)

# PARAMETERS
condition_col  <- "orig.ident"   # metadata column with condition labels
condition_1    <- "BMMC"         # first condition (ident.1)
condition_2    <- "CD34"         # second condition (ident.2)
logfc_threshold <- 0.25          # log2FC threshold for FindMarkers
test_method    <- "wilcox"       # DE test: "wilcox", "MAST", "DESeq2"
top_n_genes    <- 100            # number of top DE genes to use for enrichment
padj_cutoff    <- 0.05           # adjusted p-value cutoff for DE genes
n_pathways     <- 5              # top N pathways to display per database

# Enrichr databases to query — replace or extend as needed
databases <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Human_Gene_Atlas",
  "ENCODE_TF_ChIP-seq_2015",
  "MSigDB_Hallmark_2020"
)
# Check available databases with: enrichR::listEnrichrDbs()


# Set identities to condition
Idents(integrated) <- integrated[[condition_col, drop = TRUE]]

# Quick visualisation with built-in enrichR plot
DEenrichRPlot(
  integrated,
  ident.1         = condition_1,
  ident.2         = condition_2,
  balanced        = TRUE,
  logfc.threshold = logfc_threshold,
  test.use        = test_method,
  p.val.cutoff    = padj_cutoff,
  enrich.database = databases[1],
  num.pathway     = n_pathways,
  max.genes       = top_n_genes
) + ggtitle(paste0(databases[1], ": ", condition_1, " vs ", condition_2))

# Full DE analysis
de_results <- FindMarkers(
  integrated,
  ident.1         = condition_1,
  ident.2         = condition_2,
  logfc.threshold = logfc_threshold,
  test.use        = test_method
)

# Select top significant genes
top_genes <- de_results %>%
  filter(p_val_adj < padj_cutoff) %>%
  arrange(p_val_adj) %>%
  head(top_n_genes) %>%
  rownames()

cat("Genes submitted for enrichment:", length(top_genes), "\n")

# Run enrichment across all databases
enrichment_results <- enrichr(top_genes, databases)

# Print top pathways per database
for (db in databases) {
  cat("\nTop pathways —", db, ":\n")
  top_db <- enrichment_results[[db]] %>%
    arrange(Adjusted.P.value) %>%
    head(n_pathways)
  print(top_db[, c("Term", "Overlap", "P.value", "Adjusted.P.value")])
}

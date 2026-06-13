library(Seurat)
library(ggplot2)

# =============================================================================
# PARAMETERS
# =============================================================================
regress_cell_cycle <- FALSE   # set TRUE to regress out cell cycle effects
cc_gene_source     <- "seurat" # "seurat" uses built-in cc.genes; or provide your own

# Only used if cc_gene_source != "seurat":
s_genes_path   <- "path/to/s_genes.txt"    # one gene per line
g2m_genes_path <- "path/to/g2m_genes.txt"  # one gene per line
# =============================================================================

# Load cell cycle gene sets
if (cc_gene_source == "seurat") {
  s_genes   <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
} else {
  s_genes   <- readLines(s_genes_path)
  g2m_genes <- readLines(g2m_genes_path)
}

# Score cells for cell cycle phase
seu_filtered <- CellCycleScoring(
  seu_filtered,
  s.features   = s_genes,
  g2m.features = g2m_genes,
  set.ident    = TRUE
)

cat("Cell cycle phase distribution:\n")
print(table(seu_filtered$Phase))

# Visualise cell cycle marker distributions
RidgePlot(seu_filtered,
          features = c("PCNA", "TOP2A", "MCM6", "MKI67"),
          ncol = 2)

# PCA on cell cycle genes to check separation by phase
seu_filtered <- RunPCA(seu_filtered, features = c(s_genes, g2m_genes),
                       verbose = FALSE)
DimPlot(seu_filtered) + ggtitle("PCA coloured by cell cycle phase")

# Optional: regress out cell cycle effects during scaling
if (regress_cell_cycle) {
  cat("Regressing out cell cycle scores...\n")
  seu_filtered <- ScaleData(seu_filtered,
                            vars.to.regress = c("S.Score", "G2M.Score"),
                            features        = rownames(seu_filtered))
  # Re-run PCA to confirm cell cycle is removed
  seu_filtered <- RunPCA(seu_filtered,
                         features = VariableFeatures(seu_filtered),
                         verbose  = FALSE)
  seu_filtered <- RunPCA(seu_filtered,
                         features = c(s_genes, g2m_genes),
                         verbose  = FALSE)
  DimPlot(seu_filtered) + ggtitle("PCA after cell cycle regression")
}

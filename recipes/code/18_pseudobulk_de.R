library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)

# PARAMETERS
cell_type_col  <- "cell"       # metadata column with cell type labels
sample_col     <- "samples"    # metadata column with sample IDs
condition_col  <- "condition"  # metadata column with condition labels
assay_name     <- "RNA"        # assay to aggregate (use "RNA" not "integrated")
target_cell_type <- "B cells"  # cell type to run DE on — replace with your cell type
condition_ref  <- "Control"    # reference condition level
condition_test <- "Stimulated" # condition to compare against reference
min_counts     <- 10           # minimum total counts per gene to keep
output_file    <- "pseudobulk_DE_results.csv"


# Create sample-level labels combining cell type and sample
integrated$samples <- paste0(integrated[[condition_col, drop=TRUE]],
                              integrated[[sample_col, drop=TRUE]])

DefaultAssay(integrated) <- assay_name

# Join Seurat v5 layers if needed
if (inherits(integrated[[assay_name]], "Assay5")) {
  integrated[[assay_name]] <- JoinLayers(integrated[[assay_name]])
}

# Aggregate counts to pseudobulk level (cell type × sample)
cts <- AggregateExpression(
  integrated,
  group.by    = c(cell_type_col, "samples"),
  assays      = assay_name,
  slot        = "counts",
  return.seurat = FALSE
)[[assay_name]]

# Transpose and split by cell type
cts_t  <- as.data.frame(t(cts))
split_by <- gsub("_.*", "", rownames(cts_t))
cts_split <- split.data.frame(cts_t, f = factor(split_by))

# Clean up row names (remove cell type prefix)
cts_split <- lapply(cts_split, function(x) {
  rownames(x) <- gsub(".*_(.*)", "\\1", rownames(x))
  t(x)
})

# Extract the target cell type
counts_target <- cts_split[[target_cell_type]]
cat("Running pseudobulk DE for cell type:", target_cell_type, "\n")
cat("Samples:", ncol(counts_target), "| Genes:", nrow(counts_target), "\n")

# Build sample metadata
col_data <- data.frame(samples = colnames(counts_target)) %>%
  mutate(condition = ifelse(grepl(condition_test, samples, ignore.case = TRUE),
                            condition_test, condition_ref)) %>%
  column_to_rownames(var = "samples")
col_data$condition <- factor(col_data$condition, levels = c(condition_ref, condition_test))

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_target,
  colData   = col_data,
  design    = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= min_counts, ]
dds <- DESeq(dds)

cat("Available contrasts:", paste(resultsNames(dds), collapse = ", "), "\n")
res <- results(dds, name = paste0("condition_", condition_test, "_vs_", condition_ref))
res_df <- as.data.frame(res) %>% arrange(padj)

cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
write.csv(res_df, file = output_file, row.names = TRUE)
cat("Results saved to:", output_file, "\n")
print(head(res_df, 10))

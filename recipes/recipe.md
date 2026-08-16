---
id: scrna-seurat
name: Single-cell RNA-seq analysis with Seurat
version: 4
type: user
---


# Single-cell RNA-seq analysis with Seurat

## What this recipe does
This recipe performs single-cell RNA-seq downstream analysis in R using Seurat. It covers the full pipeline from raw count import through quality control, normalisation, dimensionality reduction, clustering, cell type annotation, and differential expression. It supports both single-sample and multi-sample workflows. The multi-sample workflow adds batch correction and integration, pseudobulk differential expression, trajectory inference, cell-cell communication, and pathway enrichment analysis.

## When to use it
- You have scRNA-seq data in any of these formats: 10X Cell Ranger output (MTX), HDF5, RDS count matrix, SingleCellExperiment object, or tab-delimited expression matrix
- You want to perform QC, normalisation, clustering, or marker gene analysis in Seurat
- You have multiple samples or conditions and want to integrate them and correct for batch effects
- You want to compare cell type compositions or gene expression across conditions
- You are asked any question about scRNA-seq analysis, clustering, UMAP, PCA, differential expression, marker genes, trajectory, or cell-cell communication

Do not use this recipe when:
- You only have raw FASTQ files (run Cell Ranger or STARsolo first)
- You need bulk RNA-seq analysis
- You want a Python-based workflow (use Scanpy instead)

## Critical execution rules
Follow these rules in every generated script unless the user explicitly states otherwise.

1. **Do not invent inputs.** Generate code only for inputs the user actually described. Do not assume an existing Seurat `.rds` object unless one is stated. Do not invent sample directories (`sample1/`, `sample2/`), metadata tables, or HDF5 files that were not mentioned.
2. **Use the working directory when that's what's available.** If the prompt names `matrix.mtx`, `genes.tsv`/`features.tsv`, and `barcodes.tsv` with no other path given, read them from the current working directory (`"."`), not a placeholder path like `"path/to/data"` or `"sample1/"`.
3. **Feature names must be unique before creating a Seurat object.** This is the single most common execution failure in this domain. Prefer `Read10X()`, which already handles this correctly — do not manually reassign feature names afterward:
   ```r
   library(Seurat)
   counts <- Read10X(data.dir = ".")
   seu <- CreateSeuratObject(counts = counts, project = "project", min.cells = 3, min.features = 200)
   ```
   Only when MTX files must be loaded manually (e.g. `Matrix::readMM()` directly), always wrap gene symbols in `make.unique()` — never assign raw, possibly-duplicated symbols as row names:
   ```r
   genes    <- read.delim("genes.tsv", header = FALSE, stringsAsFactors = FALSE)
   barcodes <- read.delim("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
   counts   <- Matrix::readMM("matrix.mtx")
   rownames(counts) <- make.unique(as.character(genes[[1]]))  # prevents "Duplicate rownames not allowed"
   colnames(counts) <- as.character(barcodes[[1]])
   ```
4. **Do not install packages inside the generated script.** Assume required packages are already available. Use `requireNamespace(..., quietly = TRUE)` to check optional packages and provide a safe fallback (see DoubletFinder pattern below) rather than calling `install.packages()` or `BiocManager::install()`.
5. **Respect object-state prerequisites.** Never call a function that depends on a prior step's output before that step has run — e.g. PCA must exist before `FindNeighbors()`/anything using PCA dimensions; a neighbour graph must exist before `FindClusters()`; clusters must exist before differential expression between clusters.
6. **Use Seurat v5-compatible APIs.** See "Seurat v5 compatibility" below.
7. **Explicitly `library()` every package whose functions you call directly.** Loading Seurat does not automatically attach every package it depends on to the search path — a function can be a real, correct choice and still fail with `could not find function "..."` if its package was never explicitly loaded. This is a common execution failure in plotting code specifically: `Seurat`'s plot functions (`VlnPlot()`, `FeatureScatter()`, `DimPlot()`, etc.) return `ggplot` objects, but customizing them with `ggtitle()`, `theme()`, `labs()`, or similar requires `library(ggplot2)` to be called explicitly — do not assume it is already attached. The same applies to pipe operators (`%>%`): `library(dplyr)` or `library(magrittr)` must be called explicitly if used.

## Inputs

### Single-sample
- 10X Cell Ranger MTX folder containing `barcodes.tsv`, `genes.tsv` or `features.tsv`, and `matrix.mtx`
- Or: HDF5 file from Cell Ranger (`filtered_feature_bc_matrix.h5`)
- Or: an existing Seurat object (`.rds`) if starting mid-pipeline
- Or: a SingleCellExperiment object or tab-delimited expression matrix

### Multi-sample
- Multiple RDS, MTX, or HDF5 files, one per sample
- Sample metadata table with columns for sample name, donor, condition, replicate, and any other grouping variables

## Outputs

### Single-sample
- QC plots (violin plots, scatter plots)
- Filtered, normalised, scaled Seurat object
- PCA plots and elbow plot
- UMAP embedding with cluster assignments
- Doublet-filtered Seurat object
- Marker gene tables (CSV) and visualisations
- Cell type annotated Seurat object

### Multi-sample
- All single-sample outputs per sample
- Integrated Seurat object (batch-corrected)
- Pseudobulk differential expression results (DESeq2)
- Trajectory plots with pseudotime
- Cell-cell communication network plots
- Pathway enrichment results (GO, KEGG, MSigDB)

## Parameters
All parameters are at the top of each script with commented defaults. Adjust before running. The agent should inspect the dataset and suggest appropriate ranges when values are not provided.

| Parameter | Default | Suggested range |
|-----------|---------|-----------------|
| `data_dir` | - | Must be set by user |
| `project_name` | `"my_project"` | Any string |
| `min_cells` | `3` | 2–5 |
| `min_features` | `200` | 100–500 |
| `min_feat` | `200` | 100–500 |
| `max_feat` | `2500` | 2000–6000 |
| `max_mt` | `5` | 5–25 (tissue dependent) |
| `n_variable` | `2000` | 1000–5000 |
| `n_pcs` | `20` | 10–30 (use elbow plot) |
| `resolution` | `0.5` | 0.1–1.5 |
| `doublet_rate` | `0.07` | 0.05–0.10 |
| `seed` | `42` | Any integer |

## Steps

### Single-sample workflow
1. Set seed and record session info for reproducibility
   Code: `code/01_reproducibility.R`

2. Load count data and create Seurat object (MTX, HDF5, RDS, SCE, or matrix)
   Code: `code/02_load_data.R`

3. Compute QC metrics - mitochondrial percentage, gene counts, UMI counts
   Code: `code/03_qc_metrics.R`

4. Filter low-quality cells based on QC thresholds
   Code: `code/04_filter_cells.R`

5. Normalise counts (log-normalisation or SCTransform)
   Code: `code/07_normalise.R`

6. Score cell cycle phases and optionally regress out cell cycle effects. Run this AFTER normalisation — `CellCycleScoring()` operates on normalised expression values, not raw counts.
   Code: `code/06_cell_cycle.R`

7. Identify highly variable genes
   Code: `code/08_variable_genes.R`

8. Scale data and regress technical covariates
   Code: `code/09_scale_data.R`

9. Linear dimensionality reduction - PCA and elbow plot
   Code: `code/10_pca.R`

10. Build kNN/SNN neighbour graph
    Code: `code/12_neighbours.R`

11. Cluster cells and tune resolution
    Code: `code/13_clustering.R`

12. Non-linear embedding - UMAP
    Code: `code/11_umap.R`

13. Detect and remove doublets using DoubletFinder. Run this AFTER steps 5–12, not before — `code/05_doublets.R` requires the object to already have normalisation, variable features, scaling, PCA, a neighbour graph, clusters, and UMAP computed (it checks for a `seurat_clusters` column and falls back to marking everything a singlet if that prerequisite is missing). Calling DoubletFinder earlier in the pipeline will fail.
    Code: `code/05_doublets.R`

14. Find marker genes per cluster
    Code: `code/14_marker_genes.R`

15. Annotate cell types (manual or automated with SingleR)
    Code: `code/15_cell_annotation.R`

**Canonical downstream order, once past filtering:** `NormalizeData → CellCycleScoring (optional) → FindVariableFeatures → ScaleData → RunPCA → FindNeighbors → FindClusters → RunUMAP → DoubletFinder → FindAllMarkers`. This matches what `code/05_doublets.R` actually requires as input, and keeps neighbours/clustering before UMAP rather than after.

### Multi-sample workflow (steps 1–15 above, then continue)
16. Load multiple samples with metadata
    Code: `code/16_load_multisample.R`

17. Integrate samples and correct for batch effects (Seurat CCA or Harmony)
    Code: `code/17_integration.R`

18. Pseudobulk differential expression analysis with DESeq2
    Code: `code/18_pseudobulk_de.R`

19. Trajectory inference with Monocle3
    Code: `code/19_trajectory.R`

20. Cell-cell communication analysis with CellChat
    Code: `code/20_cellchat.R`

21. Pathway enrichment analysis with enrichR
    Code: `code/21_pathway_enrichment.R`

## Workflows
- Single-sample graph: `single-sample-workflow.json`
- Multi-sample graph: `multi-sample-workflow.json`

## Mid-pipeline entry points
If the user already has a processed object, start at the matching step. The workflow JSON `input_type` and `output_type` fields define valid entry points. Examples:
- Have a raw Seurat object → start at step 3 (QC)
- Have a QC-passed/filtered object → start at step 5 (Normalise)
- Have a normalised object → start at step 7 (Variable genes) or step 9 (PCA)
- Have a clustered object with UMAP → start at step 13 (Doublets) or step 14 (Markers)
- Have a clustered object (no doublet filtering needed/already done) → start at step 14 (Markers) or step 15 (Annotation)
- Have per-sample clustered objects → start at step 17 (Integration)

## Code references
- `code/01_reproducibility.R` through `code/21_pathway_enrichment.R`

## Environment
Dependencies are listed in:
- `env/requirements.txt`

## Notes
- Use `^MT-` for human mitochondrial genes and `^mt-` for mouse. Set `mt_pattern` in `03_qc_metrics.R`.
- The DoubletFinder classification column name changes every run because it encodes parameter values. Scripts detect it dynamically with `grep()` - do not hardcode this column name.
- `set.seed()` is called at the top of every script. All scripts report what seed is in use for full reproducibility.
- Cell Ranger v2 output uses `genes.tsv`; v3 and later use `features.tsv`. `Read10X()` handles both automatically.
- Multi-sample integration (`17_integration.R`) automatically finds the common gene set across all samples, reports how many genes were removed, and warns if the intersection is unexpectedly small.
- CellChat (`20_cellchat.R`) automatically detects Seurat v5 Assay5 format and converts to v3 format with a clear message. The original RNA assay is not modified.
- Trajectory inference (`19_trajectory.R`) checks for required packages (monocle3, SeuratWrappers) at the start and gives installation instructions if they are missing.
- SingleR and celldex are Bioconductor packages. Install with `BiocManager::install(c("SingleR", "celldex"))`. See `env/r-packages.txt` for all dependencies.
- For doublet detection, follow `code/05_doublets.R` exactly: always call `run_doubletfinder_safe()`, never `doubletFinder()` directly. DoubletFinder is version-sensitive — use `reuse.pANN = NULL` for a fresh run, convert cluster labels with `as.character(seurat_object$seurat_clusters)`, and wrap the call in `tryCatch()`. If DoubletFinder fails or isn't installed, fall back to a `doublet_status <- "Singlet"` column so downstream execution doesn't fail.

## Seurat v5 compatibility
- Do not write directly to `object[["RNA"]]@data` or other slots — Seurat v5's `Assay5` objects use a different internal layer structure and this will error or silently do nothing.
- Do not assume an `Assay5` object has a `meta.features` slot.
- Use public accessors: `GetAssayData(object, layer = "data")` (not the old `slot = "data"` argument) to read normalised expression.
- Use `JoinLayers()` before any downstream step that needs a single combined matrix, if the assay has multiple layers (e.g. after `SCTransform` or after per-sample processing).
- Add cell-level metadata with `AddMetaData()` or `object$column <- value`, not direct slot assignment.

## Version-sensitive function rules
- **`VariableFeaturePlot()`** takes no `nfeatures`, `topfeatures`, or `label` argument. To label top genes, use `LabelPoints()` separately:
  ```r
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(seu), 10)
  p <- VariableFeaturePlot(seu)
  LabelPoints(plot = p, points = top10, repel = TRUE)
  ```
- **`CellCycleScoring()`** requires both `s.features` and `g2m.features` — never call it with only the object:
  ```r
  seu <- CellCycleScoring(seu, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)
  ```
- **DoubletFinder**: use the single, currently-installed API convention exposed by the package (see `code/05_doublets.R`'s `run_doubletfinder_safe()`). Do not mix `doubletFinder_v3()`-style and `doubletFinder()`-style calling conventions, and do not invent parameters not shown in that script.

## Pre-output self-check
Before returning generated code, verify:
- Every referenced input file was actually supplied or explicitly stated as available — nothing invented.
- No placeholder path (`path/to/data`, `sample1/`) remains when a working-directory input is available.
- Feature/gene names are unique before `CreateSeuratObject()`.
- Required packages are loaded; none are installed inside the script.
- PCA exists before anything using PCA dimensions; a neighbour graph exists before clustering; clusters exist before differential expression.
- The final object is assigned to a clearly discoverable variable.

## Limitations
- **Ambient RNA removal is not included.** SoupX and DecontX require the raw unfiltered Cell Ranger output alongside the filtered matrix. Not all users have this, and the correction is dataset-specific. This is a genuine methodological gap - consider running SoupX before this recipe if ambient RNA contamination is a concern for your tissue type.
- **Single modality only.** This recipe handles gene expression (RNA) data only. ATAC-seq, protein (CITE-seq), or spatial transcriptomics data require different workflows and are not covered here.
- **No automated parameter optimisation.** QC thresholds, number of PCs, and clustering resolution must be set manually. The recipe provides recommended ranges and plots to guide these choices, but optimal values depend on your specific dataset and biology.
- **CellChat and trajectory inference are computationally intensive.** Datasets with more than 50,000 cells may require substantial RAM (32GB+) and run time. Consider downsampling or running on an HPC for large datasets.
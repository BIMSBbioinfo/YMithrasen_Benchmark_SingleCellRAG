# Installs the R packages the benchmark needs and prints their versions.
#   Rscript setup_env.R

# DoubletFinder 2.0.6
DOUBLETFINDER_REF <- "1b244d8f0d54b4b1cb4365639931bbb16f01e1cd"

# compile packages on all cores but one
options(Ncpus = max(1L, parallel::detectCores() - 1L))

cran <- c("Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "Matrix",
          "tidyr", "stringr", "readr", "ggrepel", "jsonlite", "remotes", "BiocManager")
bioc <- c("SingleCellExperiment", "scDblFinder", "SingleR", "celldex")

if (is.na(getOption("repos")["CRAN"]) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

not_installed <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

# versions the benchmark was run with; installed first so the rest builds against them
pinned <- c(Matrix = "1.7-4", readr = "2.1.6", SeuratObject = "5.4.0", Seurat = "5.5.0")

if (length(not_installed("remotes"))) install.packages("remotes")
for (p in names(pinned)) {
  have <- if (requireNamespace(p, quietly = TRUE)) packageVersion(p) else NULL
  if (is.null(have) || have != pinned[[p]]) {
    tryCatch(
      remotes::install_version(p, version = pinned[[p]], upgrade = "never"),
      error = function(e) message("could not install ", p, " ", pinned[[p]], ": ", conditionMessage(e))
    )
  }
}

if (length(not_installed(cran))) install.packages(not_installed(cran))
if (length(not_installed(bioc))) BiocManager::install(not_installed(bioc), update = FALSE, ask = FALSE)

# DoubletFinder is only on GitHub
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github(
    paste0("chris-mcginnis-ucsf/DoubletFinder@", DOUBLETFINDER_REF),
    upgrade = "never"
  )
}

required <- setdiff(c(cran, bioc, "DoubletFinder"), c("remotes", "BiocManager"))
missing <- not_installed(required)
if (length(missing)) {
  cat("MISSING PACKAGES:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}

cat(R.version.string, "\n")
for (p in required) cat(sprintf("%-22s %s\n", p, as.character(packageVersion(p))))
cat("All required R packages are installed.\n")

for (p in names(pinned)) {
  if (packageVersion(p) != pinned[[p]]) {
    cat(sprintf("NOTE: %s is %s, the benchmark was run with %s\n", p, as.character(packageVersion(p)), pinned[[p]]))
  }
}
#!/usr/bin/env Rscript
# =============================================================================
# install_packages.R -- install the R packages used by the pipeline.
# Only missing packages are installed. Usage: Rscript install_packages.R
# =============================================================================

# Use the CRAN mirror already configured (e.g. Posit Package Manager on CI);
# fall back to the cloud mirror.
cran_repo <- unname(getOption("repos")["CRAN"])
if (length(cran_repo) == 0 || is.na(cran_repo) || cran_repo == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
ncpus <- max(1L, parallel::detectCores() - 1L)

cran <- c(
  "yaml", "data.table", "R.utils", "tidyverse", "matrixStats", "pheatmap",
  "RColorBrewer", "ggrepel", "scales", "patchwork", "viridis", "hexbin",
  "msigdbr", "rmarkdown", "knitr", "remotes", "BiocManager",
  # xCell dependencies (installed from GitHub below)
  "pracma", "quadprog", "digest", "curl", "MASS"
)
bioc <- c(
  "DESeq2", "apeglm", "GEOquery", "Biobase", "SummarizedExperiment",
  "AnnotationDbi", "org.Hs.eg.db", "GO.db", "clusterProfiler", "DOSE",
  "enrichplot", "ReactomePA", "reactome.db", "GSVA", "GSEABase"
)

missing <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install <- missing(c(cran, bioc))
if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  BiocManager::install(to_install, update = FALSE, ask = FALSE, Ncpus = ncpus)
}

if (!requireNamespace("xCell", quietly = TRUE)) {
  remotes::install_github("dviraran/xCell", dependencies = FALSE, upgrade = "never")
}

still_missing <- missing(c(cran, bioc, "xCell"))
if (length(still_missing) > 0) {
  stop("Packages failed to install: ", paste(still_missing, collapse = ", "))
}
message("All packages available. R ", getRversion(),
        ", Bioconductor ", as.character(BiocManager::version()))

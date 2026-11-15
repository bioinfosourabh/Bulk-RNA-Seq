#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R -- run the numbered scripts in order.
#
# Usage (from the repository root):
#   Rscript run_pipeline.R                              # config.yaml, steps 00A-07
#   Rscript run_pipeline.R --download                   # also run 00_download_data.R first
#   Rscript run_pipeline.R --config example/config.yaml # simulated example data
#   Rscript run_pipeline.R --from 02 --to 07            # a subset of steps
#
# Each step runs in a fresh R process; the run stops at the first failure.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) > i) return(args[i + 1])
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  default
}

config <- get_arg("--config", "config.yaml")
from <- get_arg("--from", NULL)
to <- get_arg("--to", NULL)

steps <- c(
  "00_download_data.R",
  "00A_make_counts_from_processed_matrix.R",
  "00B_make_coldata_from_GEO.R",
  "01_input_validation.R",
  "02_deseq2_analysis.R",
  "03_qc_vst_pca_distances.R",
  "04_deconvolution_xcell.R",
  "05_enrichment_ORA_GSEA.R",
  "07_visualizations.R",
  "06_render_report.R"   # renders the report after the figures exist
)
if (!"--download" %in% args) steps <- steps[-1]
ids <- sub("_.*$", "", steps)
if (!is.null(from)) steps <- steps[seq(match(from, ids), length(steps))]
ids <- sub("_.*$", "", steps)
if (!is.null(to)) steps <- steps[seq_len(match(to, ids))]
if (!file.exists(config)) stop("Configuration not found: ", config)

rscript <- file.path(R.home("bin"), "Rscript")
t0 <- Sys.time()
for (s in steps) {
  cat("\n>>> ", s, "\n", sep = "")
  t1 <- Sys.time()
  status <- system2(rscript, c(file.path("scripts", s), "--config", config))
  if (status != 0) stop("Step failed: ", s, " (exit status ", status, ")")
  cat(sprintf(">>> %s finished in %.1f s\n", s, as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}
cat(sprintf("\nPipeline finished in %.1f min (%d steps, config: %s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(steps), config))

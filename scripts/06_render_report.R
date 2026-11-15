#!/usr/bin/env Rscript
# =============================================================================
# 06_render_report.R
#
# Render reports/report.Rmd into a self-contained HTML summary of the run
# (sample counts, filtering, QC flags, DE summary and top genes, xCell and
# enrichment summaries, key figures). Requires pandoc (bundled with RStudio;
# otherwise install it separately).
#
# Output: results/report.html
# =============================================================================

source("scripts/utils.R")

cfg <- load_config()
log_step("06. Render HTML report")

if (!rmarkdown::pandoc_available()) {
  stop("pandoc is not available; install pandoc to render the report.")
}
require_files(p_table(cfg, "DESeq2_results_all.tsv"))

results_dir <- normalizePath(cfg$dirs$results)
out <- rmarkdown::render(
  input = "reports/report.Rmd",
  output_file = "report.html",
  output_dir = results_dir,
  intermediates_dir = tempdir(),
  knit_root_dir = getwd(),
  params = list(title = cfg$title, results_dir = results_dir,
                config_file = basename(cfg$config_file),
                reference = cfg$design$reference, treatment = cfg$design$treatment,
                alpha = cfg$deseq2$alpha),
  envir = new.env(),
  quiet = TRUE
)
log_msg("Report written to ", out)

write_session_info(cfg, "06_render_report")

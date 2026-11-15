#!/usr/bin/env Rscript
# =============================================================================
# 00_download_data.R
#
# Download GSE192341 from NCBI GEO:
#   * the series matrix (sample annotation), exported to data_raw/<acc>_pdata.tsv
#   * the supplementary processed count matrix, saved as the file named in
#     config$input$processed_matrix
#
# Requires internet access to ftp.ncbi.nlm.nih.gov. Not needed for the example
# run (example/config.yaml), which ships its own inputs.
#
# Usage: Rscript scripts/00_download_data.R [--config config.yaml]
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

cfg <- load_config()
acc <- cfg$geo$accession
log_step(paste("00. Download", acc, "from GEO"))

options(timeout = max(1200, getOption("timeout")))
dl_dir <- file.path(cfg$dirs$data_raw, "geo_download")
dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Sample annotation from the series matrix ------------------------------
log_msg("Fetching series matrix")
gse <- getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE, destdir = dl_dir)
if (length(gse) > 1) {
  log_msg("Series has ", length(gse), " platforms; combining sample annotation")
}
pdata <- do.call(rbind, lapply(gse, function(es) {
  pd <- pData(es)
  pd[] <- lapply(pd, as.character)
  pd
}))
pdata_file <- cfg$input$pdata
dir.create(dirname(pdata_file), recursive = TRUE, showWarnings = FALSE)
write_tsv_base(pdata, pdata_file)
log_msg("Sample annotation: ", nrow(pdata), " samples x ", ncol(pdata),
        " fields -> ", pdata_file)

# ---- Supplementary processed matrix ----------------------------------------
log_msg("Fetching supplementary files")
supp <- getGEOSuppFiles(acc, makeDirectory = FALSE, baseDir = dl_dir)
supp_files <- rownames(supp)
log_msg("Downloaded:\n  ", paste(basename(supp_files), collapse = "\n  "))

hit <- supp_files[grepl(cfg$geo$processed_file_regex, basename(supp_files),
                        ignore.case = TRUE)]
hit <- hit[!grepl("\\.tar$", hit)]
if (length(hit) == 0) {
  stop("No supplementary file matched '", cfg$geo$processed_file_regex,
       "'. Inspect ", dl_dir, " and set geo$processed_file_regex or ",
       "input$processed_matrix in the configuration.")
}
src <- hit[1]
dest <- cfg$input$processed_matrix

# Decompress if needed so that 00A can read the file directly.
if (grepl("\\.gz$", src) && !grepl("\\.gz$", dest)) {
  log_msg("Decompressing ", basename(src))
  con_in <- gzfile(src, open = "rt")
  con_out <- file(dest, open = "wt")
  repeat {
    lines <- readLines(con_in, n = 10000)
    if (length(lines) == 0) break
    writeLines(lines, con_out)
  }
  close(con_in)
  close(con_out)
} else {
  file.copy(src, dest, overwrite = TRUE)
}
log_msg("Processed matrix -> ", dest)

# ---- Manifest ----------------------------------------------------------------
manifest <- data.frame(
  file = basename(supp_files),
  size_bytes = file.size(supp_files),
  md5 = unname(tools::md5sum(supp_files)),
  downloaded = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
write_tsv_base(manifest, file.path(cfg$dirs$data_raw, "download_manifest.tsv"))

write_session_info(cfg, "00_download_data")

#!/usr/bin/env Rscript
# =============================================================================
# 00B_make_coldata_from_GEO.R
#
# Build DESeq2 sample annotation (coldata) from the GEO series annotation
# exported by 00_download_data.R.
#
# Input : config$input$pdata (one row per GEO sample)
# Output: data/coldata.tsv, results/rds/coldata.rds
# =============================================================================

source("scripts/utils.R")

cfg <- load_config()
log_step("00B. Sample metadata construction")

require_files(cfg$input$pdata)
pdata <- read_tsv_base(cfg$input$pdata)
log_msg("GEO annotation: ", nrow(pdata), " samples")

ch1_cols <- grep(":ch1$", colnames(pdata), value = TRUE)

# ---- Response column ----------------------------------------------------------
resp_col <- cfg$metadata$response_column
if (is.null(resp_col)) {
  cand <- ch1_cols[grepl("pcr|response", ch1_cols, ignore.case = TRUE)]
  if (length(cand) == 0) {
    stop("Could not detect the response column. Available ':ch1' columns:\n  ",
         paste(ch1_cols, collapse = "\n  "),
         "\nSet metadata$response_column in the configuration.")
  }
  resp_col <- cand[1]
}
if (!resp_col %in% colnames(pdata)) stop("Column not found: ", resp_col)
log_msg("Response column: '", resp_col, "'")

raw_resp <- trimws(tolower(as.character(pdata[[resp_col]])))
ref <- cfg$design$reference
trt <- cfg$design$treatment
condition <- ifelse(raw_resp %in% tolower(cfg$metadata$pcr_values), trt,
             ifelse(raw_resp %in% tolower(cfg$metadata$no_pcr_values), ref, NA))

cat("\nMapping of GEO values to condition:\n")
print(table(GEO = pdata[[resp_col]], condition = condition, useNA = "ifany"))

# ---- Assemble coldata ----------------------------------------------------------
coldata <- data.frame(
  sample_id = pdata$geo_accession,
  title = pdata$title,
  response = pdata[[resp_col]],
  stringsAsFactors = FALSE
)
for (nm in names(cfg$metadata$extra_columns)) {
  pattern <- cfg$metadata$extra_columns[[nm]]
  hit <- ch1_cols[grepl(pattern, ch1_cols, ignore.case = TRUE) & ch1_cols != resp_col]
  if (length(hit) > 0) {
    coldata[[nm]] <- pdata[[hit[1]]]
    log_msg("Added ", nm, " from '", hit[1], "'")
  }
}
coldata[[cfg$design$condition_column]] <- factor(condition, levels = c(ref, trt))
rownames(coldata) <- coldata$sample_id

cat("\nSamples per condition:\n")
print(table(coldata[[cfg$design$condition_column]], useNA = "ifany"))

write_tsv_base(coldata, p_data(cfg, "coldata.tsv"))
saveRDS(coldata, p_rds(cfg, "coldata.rds"))
log_msg("Wrote ", p_data(cfg, "coldata.tsv"))

write_session_info(cfg, "00B_make_coldata_from_GEO")

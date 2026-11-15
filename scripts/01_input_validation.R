#!/usr/bin/env Rscript
# =============================================================================
# 01_input_validation.R
#
# Align count-matrix columns with the sample annotation, drop samples without
# a response label, compute per-sample QC metrics and remove low-count genes.
#
# Input : results/rds/counts_raw.rds, results/rds/coldata.rds
# Output: results/rds/counts_filtered.rds, results/rds/coldata_validated.rds,
#         results/tables/sample_qc_metrics.tsv,
#         results/tables/gene_filtering_summary.tsv
# =============================================================================

source("scripts/utils.R")

cfg <- load_config()
log_step("01. Input validation and gene filtering")

require_files(p_rds(cfg, "counts_raw.rds"), p_rds(cfg, "coldata.rds"))
counts <- readRDS(p_rds(cfg, "counts_raw.rds"))
coldata <- readRDS(p_rds(cfg, "coldata.rds"))
cond_col <- cfg$design$condition_column

# ---- Match count columns to coldata rows ---------------------------------------
# GEO processed matrices are often labelled with sample titles or patient codes
# rather than GSM accessions. Use the coldata column that matches the most
# count-matrix column names, and rename the columns to GSM accessions.
id_cols <- names(coldata)[vapply(coldata, function(x) is.character(x) || is.factor(x),
                                 logical(1))]
overlap <- vapply(id_cols, function(cn) sum(colnames(counts) %in% as.character(coldata[[cn]])),
                  numeric(1))
best <- names(which.max(overlap))
if (max(overlap) == 0) {
  stop("No count-matrix column names match any identifier column in coldata ",
       "(checked: ", paste(id_cols, collapse = ", "), ").")
}
log_msg("Count columns matched to coldata$", best, " (", max(overlap), " of ",
        ncol(counts), " columns)")
if (best != "sample_id") {
  key <- as.character(coldata[[best]])
  if (anyDuplicated(key[key %in% colnames(counts)])) {
    stop("coldata$", best, " is not unique; cannot map columns unambiguously.")
  }
  matched <- colnames(counts) %in% key
  colnames(counts)[matched] <- coldata$sample_id[match(colnames(counts)[matched], key)]
}

only_counts <- setdiff(colnames(counts), rownames(coldata))
only_coldata <- setdiff(rownames(coldata), colnames(counts))
if (length(only_counts) > 0) {
  warning("Samples in counts but not in coldata (dropped): ",
          paste(only_counts, collapse = ", "))
}
if (length(only_coldata) > 0) {
  warning("Samples in coldata but not in counts (dropped): ",
          paste(only_coldata, collapse = ", "))
}

common <- intersect(colnames(counts), rownames(coldata))
counts <- counts[, common, drop = FALSE]
coldata <- coldata[common, , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(coldata)))

# ---- Remove samples without a condition label -----------------------------------
has_label <- !is.na(coldata[[cond_col]])
if (any(!has_label)) {
  log_msg("Removing ", sum(!has_label), " sample(s) with missing ", cond_col, ": ",
          paste(rownames(coldata)[!has_label], collapse = ", "))
}
counts <- counts[, has_label, drop = FALSE]
coldata <- coldata[has_label, , drop = FALSE]
coldata[[cond_col]] <- droplevels(coldata[[cond_col]])

group_sizes <- table(coldata[[cond_col]])
print(group_sizes)
if (length(group_sizes) != 2 || any(group_sizes < 2)) {
  stop("Two conditions with at least two samples each are required.")
}

# ---- Per-sample QC metrics -----------------------------------------------------
sample_qc <- data.frame(
  sample = colnames(counts),
  condition = as.character(coldata[[cond_col]]),
  library_size = colSums(counts),
  detected_genes = colSums(counts > 0),
  pct_top100_genes = apply(counts, 2, function(x) {
    100 * sum(sort(x, decreasing = TRUE)[1:min(100, length(x))]) / sum(x)
  })
)
write_tsv_base(sample_qc, p_table(cfg, "sample_qc_metrics.tsv"))

# ---- Gene filtering --------------------------------------------------------------
min_count <- cfg$filtering$min_count
min_samples <- cfg$filtering$min_samples
if (identical(min_samples, "smallest_group")) min_samples <- min(group_sizes)
min_samples <- as.integer(min_samples)

keep <- rowSums(counts >= min_count) >= min_samples
counts_filtered <- counts[keep, , drop = FALSE]
log_msg(sprintf("Gene filter: >= %d counts in >= %d samples", min_count, min_samples))
log_msg("Genes before filtering: ", nrow(counts))
log_msg("Genes after filtering:  ", nrow(counts_filtered))

write_tsv_base(data.frame(
  min_count = min_count, min_samples = min_samples,
  genes_before = nrow(counts), genes_after = nrow(counts_filtered),
  samples = ncol(counts_filtered)
), p_table(cfg, "gene_filtering_summary.tsv"))

saveRDS(counts_filtered, p_rds(cfg, "counts_filtered.rds"))
saveRDS(coldata, p_rds(cfg, "coldata_validated.rds"))

write_session_info(cfg, "01_input_validation")

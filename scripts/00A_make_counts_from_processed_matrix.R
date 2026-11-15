#!/usr/bin/env Rscript
# =============================================================================
# 00A_make_counts_from_processed_matrix.R
#
# Parse the GEO processed matrix into a gene x sample count matrix for DESeq2.
#
# Input : config$input$processed_matrix
#         (col 1 gene id, optional non-numeric annotation columns, then one
#          numeric column per sample)
# Output: data/counts.csv, data/gene_annotation.csv, results/rds/counts_raw.rds
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(data.table)
})

cfg <- load_config()
log_step("00A. Count matrix construction")

infile <- cfg$input$processed_matrix
require_files(infile)
raw <- data.table::fread(infile, data.table = FALSE, check.names = FALSE)
log_msg("Read ", nrow(raw), " rows x ", ncol(raw), " columns from ", infile)

# ---- Identify annotation and sample columns ---------------------------------
is_num <- vapply(raw, is.numeric, logical(1))
if (is_num[1]) stop("The first column must contain gene identifiers.")
sample_cols <- colnames(raw)[is_num]
annot_cols <- colnames(raw)[!is_num]
if (length(sample_cols) < 2) stop("Fewer than two numeric sample columns found.")
log_msg("Annotation columns: ", paste(annot_cols, collapse = ", "))
log_msg("Sample columns: ", length(sample_cols))

gene_id <- strip_ensembl_version(as.character(raw[[1]]))
gene_symbol <- if (length(annot_cols) >= 2) as.character(raw[[annot_cols[2]]]) else NA_character_

counts <- as.matrix(raw[, sample_cols, drop = FALSE])
rownames(counts) <- gene_id

# ---- Sanity checks ------------------------------------------------------------
if (anyNA(counts)) stop("The matrix contains missing values.")
if (any(counts < 0)) stop("The matrix contains negative values; raw counts are required.")

non_int <- mean(abs(counts - round(counts)) > 1e-8)
if (non_int > 0) {
  med_lib <- stats::median(colSums(counts))
  if (non_int > 0.5 && med_lib < 2e6) {
    stop("Values look like normalised expression (non-integer, median column ",
         "sum ", round(med_lib), "). DESeq2 requires raw or expected counts.")
  }
  log_msg(sprintf("%.1f%% of values are non-integer (e.g. expected counts); rounding.",
                  100 * non_int))
  counts <- round(counts)
}
storage.mode(counts) <- "integer"

# ---- Duplicated identifiers (e.g. PAR_Y genes after version stripping) --------
if (anyDuplicated(rownames(counts))) {
  n_dup <- sum(duplicated(rownames(counts)))
  log_msg("Summing ", n_dup, " rows with duplicated gene identifiers")
  counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

gene_annotation <- unique(data.frame(gene_id = gene_id, gene_symbol = gene_symbol,
                                     stringsAsFactors = FALSE))
gene_annotation <- gene_annotation[!duplicated(gene_annotation$gene_id), ]

log_msg("Count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")
log_msg("Library sizes (millions): ",
        paste(range(round(colSums(counts) / 1e6, 2)), collapse = " - "))

utils::write.csv(counts, p_data(cfg, "counts.csv"), row.names = TRUE)
utils::write.csv(gene_annotation, p_data(cfg, "gene_annotation.csv"), row.names = FALSE)
saveRDS(counts, p_rds(cfg, "counts_raw.rds"))

write_session_info(cfg, "00A_make_counts_from_processed_matrix")

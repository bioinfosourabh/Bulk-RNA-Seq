#!/usr/bin/env Rscript
# =============================================================================
# 04_deconvolution_xcell.R
#
# xCell enrichment scores for 64 immune and stromal cell types (plus the
# ImmuneScore, StromaScore and MicroenvironmentScore composites), followed by
# a two-sided Wilcoxon rank-sum test per cell type with BH correction.
#
# xCell is rank-based within each sample, so the monotone VST transform gives
# the same within-sample ranks as size-factor-normalised counts. Counts are not
# corrected for gene length (xCell was trained on TPM/FPKM-like data); scores
# should be read as relative enrichment between samples, not as proportions.
#
# Input : results/rds/vsd.rds, results/rds/coldata_validated.rds
# Output: results/tables/deconvolution_xcell_scores.tsv
#         results/tables/deconvolution_xcell_pcr_vs_nopcr.tsv
#         results/rds/xcell_scores.rds
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(xCell)
})

cfg <- load_config()
log_step("04. Cellular deconvolution (xCell)")

if (!isTRUE(cfg$xcell$run)) {
  log_msg("xcell$run is false; skipping.")
  quit(save = "no", status = 0)
}

require_files(p_rds(cfg, "vsd.rds"), p_rds(cfg, "coldata_validated.rds"))
vsd <- readRDS(p_rds(cfg, "vsd.rds"))
coldata <- readRDS(p_rds(cfg, "coldata_validated.rds"))
cond_col <- cfg$design$condition_column
ref <- cfg$design$reference
trt <- cfg$design$treatment

# ---- Gene symbols as row names, one row per symbol -----------------------------
expr <- assay(vsd)
symbols <- map_to_symbol(rownames(expr), cfg)
ok <- !is.na(symbols) & symbols != "" & !grepl("^ENS", symbols)
expr <- expr[ok, , drop = FALSE]
symbols <- symbols[ok]
# For duplicated symbols keep the row with the highest mean expression.
ord <- order(rowMeans(expr), decreasing = TRUE)
expr <- expr[ord, , drop = FALSE]
symbols <- symbols[ord]
keep <- !duplicated(symbols)
expr <- expr[keep, , drop = FALSE]
rownames(expr) <- symbols[keep]

data("xCell.data", package = "xCell", envir = environment())
n_shared <- length(intersect(rownames(expr), xCell.data$genes))
log_msg("Genes with symbols: ", nrow(expr), "; shared with xCell signatures: ",
        n_shared, " of ", length(xCell.data$genes))
if (n_shared < 1000) {
  warning("Fewer than 1,000 xCell signature genes are present; scores may be unreliable.")
}

# ---- xCell -------------------------------------------------------------------
xcell_scores <- xCellAnalysis(expr, rnaseq = TRUE, parallel.sz = cfg$xcell$parallel_sz)
log_msg("xCell scores: ", nrow(xcell_scores), " cell types x ", ncol(xcell_scores), " samples")

xcell_df <- data.frame(sample = colnames(xcell_scores), t(xcell_scores),
                       check.names = FALSE)
xcell_df$condition <- as.character(coldata[xcell_df$sample, cond_col])
xcell_df <- xcell_df[, c("sample", "condition",
                         setdiff(names(xcell_df), c("sample", "condition")))]

# ---- Group comparison ----------------------------------------------------------
is_trt <- xcell_df$condition == trt
is_ref <- xcell_df$condition == ref
comparison <- do.call(rbind, lapply(rownames(xcell_scores), function(ct) {
  x <- xcell_scores[ct, xcell_df$sample]
  p <- if (stats::sd(x) > 0) {
    suppressWarnings(stats::wilcox.test(x[is_trt], x[is_ref], exact = FALSE)$p.value)
  } else NA_real_
  data.frame(population = ct,
             mean_treatment = mean(x[is_trt]),
             mean_reference = mean(x[is_ref]),
             median_treatment = stats::median(x[is_trt]),
             median_reference = stats::median(x[is_ref]),
             difference = mean(x[is_trt]) - mean(x[is_ref]),
             p_value = p)
}))
comparison$padj <- stats::p.adjust(comparison$p_value, method = "BH")
comparison <- comparison[order(comparison$padj, comparison$p_value, na.last = TRUE), ]
names(comparison) <- sub("treatment", trt, sub("reference", ref, names(comparison)))

log_msg("Cell types with padj < 0.05: ", sum(comparison$padj < 0.05, na.rm = TRUE))

write_tsv_base(xcell_df, p_table(cfg, "deconvolution_xcell_scores.tsv"))
write_tsv_base(comparison, p_table(cfg, "deconvolution_xcell_pcr_vs_nopcr.tsv"))
saveRDS(xcell_scores, p_rds(cfg, "xcell_scores.rds"))

write_session_info(cfg, "04_deconvolution_xcell")

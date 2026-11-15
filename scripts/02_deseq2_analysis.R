#!/usr/bin/env Rscript
# =============================================================================
# 02_deseq2_analysis.R
#
# Differential expression (treatment vs reference) with DESeq2, apeglm
# log2 fold-change shrinkage and Benjamini-Hochberg FDR.
#
# Input : results/rds/counts_filtered.rds, results/rds/coldata_validated.rds
# Output: results/rds/dds.rds
#         results/tables/DESeq2_results_all.tsv
#         results/tables/DESeq2_results_significant.tsv
#         results/tables/normalized_counts.tsv, results/tables/vst_matrix.tsv
#         results/tables/DESeq2_run_summary.tsv
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

cfg <- load_config()
log_step("02. Differential expression analysis (DESeq2 + apeglm)")

require_files(p_rds(cfg, "counts_filtered.rds"), p_rds(cfg, "coldata_validated.rds"))
counts <- readRDS(p_rds(cfg, "counts_filtered.rds"))
coldata <- readRDS(p_rds(cfg, "coldata_validated.rds"))

cond_col <- cfg$design$condition_column
ref <- cfg$design$reference
trt <- cfg$design$treatment
alpha <- cfg$deseq2$alpha

coldata[[cond_col]] <- stats::relevel(factor(coldata[[cond_col]]), ref = ref)
design <- stats::as.formula(cfg$design$formula)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design)
dds <- DESeq(dds)

coef_name <- paste0(cond_col, "_", trt, "_vs_", ref)
log_msg("Model coefficients: ", paste(resultsNames(dds), collapse = ", "))
if (!coef_name %in% resultsNames(dds)) {
  stop("Coefficient '", coef_name, "' not found. Check design$reference, ",
       "design$treatment and that the condition is the last term of the formula.")
}

# Unshrunken (MLE) results: Wald statistic, p-values and BH-adjusted p-values.
res <- results(dds, name = coef_name, alpha = alpha)
summary(res)

# apeglm shrinkage of log2 fold changes (Zhu, Ibrahim & Love 2019). The
# p-values and adjusted p-values are carried over unchanged from `res`.
res_shrunk <- lfcShrink(dds, coef = coef_name, type = cfg$deseq2$shrinkage, res = res)

res_df <- data.frame(
  gene_id = rownames(res),
  baseMean = res$baseMean,
  log2FoldChange = res_shrunk$log2FoldChange,   # apeglm-shrunken
  lfcSE = res_shrunk$lfcSE,
  log2FoldChange_MLE = res$log2FoldChange,      # unshrunken
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  stringsAsFactors = FALSE
)
res_df$symbol <- map_to_symbol(res_df$gene_id, cfg)
res_df$regulation <- ifelse(!is.na(res_df$padj) & res_df$padj < alpha,
                            ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
                            "Not Significant")
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
res_df <- res_df[, c("gene_id", "symbol", setdiff(names(res_df), c("gene_id", "symbol")))]

sig_df <- res_df[res_df$regulation != "Not Significant", ]

log_msg("Genes tested:            ", nrow(res_df))
log_msg("padj < ", alpha, ":            ", nrow(sig_df))
log_msg("  higher in ", trt, ":        ", sum(sig_df$regulation == "Upregulated"))
log_msg("  lower in ", trt, ":         ", sum(sig_df$regulation == "Downregulated"))

# ---- Export ------------------------------------------------------------------
normalized_counts <- counts(dds, normalized = TRUE)
vst_matrix <- assay(vst(dds, blind = FALSE))

saveRDS(dds, p_rds(cfg, "dds.rds"))
saveRDS(res_shrunk, p_rds(cfg, "res_shrunk.rds"))
write_tsv_base(res_df, p_table(cfg, "DESeq2_results_all.tsv"))
write_tsv_base(sig_df, p_table(cfg, "DESeq2_results_significant.tsv"))
write_tsv_base(data.frame(gene_id = rownames(normalized_counts), normalized_counts,
                          check.names = FALSE), p_table(cfg, "normalized_counts.tsv"))
write_tsv_base(data.frame(gene_id = rownames(vst_matrix), vst_matrix, check.names = FALSE),
               p_table(cfg, "vst_matrix.tsv"))
write_tsv_base(data.frame(
  contrast = paste(trt, "vs", ref), design = cfg$design$formula,
  samples = ncol(dds), genes_tested = nrow(res_df), alpha = alpha,
  significant = nrow(sig_df),
  up = sum(sig_df$regulation == "Upregulated"),
  down = sum(sig_df$regulation == "Downregulated"),
  shrinkage = cfg$deseq2$shrinkage
), p_table(cfg, "DESeq2_run_summary.tsv"))

write_session_info(cfg, "02_deseq2_analysis")

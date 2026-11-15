#!/usr/bin/env Rscript
# =============================================================================
# 03_qc_vst_pca_distances.R
#
# Sample-level QC on variance-stabilised data: PCA on the most variable genes,
# Euclidean distances, Pearson correlations and distance-based outlier flags.
# Outliers are flagged, not removed.
#
# Input : results/rds/dds.rds
# Output: results/rds/vsd.rds
#         results/tables/qc_pca_coordinates.tsv, qc_pca_variance_explained.tsv
#         results/tables/qc_sample_distance_matrix.tsv
#         results/tables/qc_sample_correlation_matrix.tsv
#         results/tables/qc_outlier_flags.tsv, qc_top_variable_genes.tsv
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
})

cfg <- load_config()
log_step("03. QC: VST, PCA, sample distances")

require_files(p_rds(cfg, "dds.rds"))
dds <- readRDS(p_rds(cfg, "dds.rds"))
cond_col <- cfg$design$condition_column

# Blind VST so that the QC is not informed by the design.
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

# ---- PCA on the top variable genes (as in DESeq2::plotPCA) -------------------
n_top <- min(cfg$qc$n_top_variable, nrow(mat))
rv <- rowVars(mat)
top_idx <- order(rv, decreasing = TRUE)[seq_len(n_top)]
pca <- prcomp(t(mat[top_idx, ]), center = TRUE, scale. = FALSE)
percent_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
n_pc <- min(10, ncol(pca$x))

pca_df <- data.frame(sample_id = colnames(vsd),
                     condition = as.character(colData(vsd)[[cond_col]]),
                     pca$x[, seq_len(min(5, ncol(pca$x))), drop = FALSE],
                     check.names = FALSE)
log_msg(sprintf("PCA on %d genes: PC1 %.1f%%, PC2 %.1f%%", n_top,
                percent_var[1], percent_var[2]))

# ---- Distances and correlations (all filtered genes) -------------------------
dist_mat <- as.matrix(dist(t(mat), method = "euclidean"))
cor_mat <- cor(mat, method = "pearson")

# ---- Outlier detection: z-score of each sample's mean distance -------------
mean_dist <- rowSums(dist_mat) / (ncol(dist_mat) - 1)
z <- as.numeric(scale(mean_dist))
outliers <- data.frame(sample_id = names(mean_dist),
                       condition = as.character(colData(vsd)[[cond_col]]),
                       mean_distance = mean_dist,
                       z_score = z,
                       outlier_flag = abs(z) > cfg$qc$outlier_z)
n_out <- sum(outliers$outlier_flag)
log_msg("Samples flagged (|z| > ", cfg$qc$outlier_z, "): ", n_out,
        if (n_out > 0) paste0(" [", paste(outliers$sample_id[outliers$outlier_flag],
                                          collapse = ", "), "]") else "")

# ---- Export ------------------------------------------------------------------
mat_to_df <- function(m) data.frame(sample_id = rownames(m), m, check.names = FALSE)

write_tsv_base(pca_df, p_table(cfg, "qc_pca_coordinates.tsv"))
write_tsv_base(data.frame(PC = paste0("PC", seq_len(n_pc)),
                          variance_explained = percent_var[seq_len(n_pc)] / 100),
               p_table(cfg, "qc_pca_variance_explained.tsv"))
write_tsv_base(mat_to_df(dist_mat), p_table(cfg, "qc_sample_distance_matrix.tsv"))
write_tsv_base(mat_to_df(cor_mat), p_table(cfg, "qc_sample_correlation_matrix.tsv"))
write_tsv_base(outliers, p_table(cfg, "qc_outlier_flags.tsv"))
write_tsv_base(data.frame(gene_id = rownames(mat)[top_idx], variance = rv[top_idx]),
               p_table(cfg, "qc_top_variable_genes.tsv"))
saveRDS(vsd, p_rds(cfg, "vsd.rds"))

write_session_info(cfg, "03_qc_vst_pca_distances")

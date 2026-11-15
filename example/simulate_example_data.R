#!/usr/bin/env Rscript
# =============================================================================
# simulate_example_data.R
#
# Generates the small simulated dataset in example/data_raw/ that is used to
# test the pipeline end-to-end (example/config.yaml, GitHub Actions CI).
# The files are committed, so this script only needs to be re-run to change
# the simulation. It is NOT real patient data.
#
# Design: 18 "tumours", 7 pCR, 10 No pCR and 1 with missing response (to
# exercise sample exclusion). Counts are negative binomial. Three correlated
# gene modules vary between samples:
#   * immune / interferon-gamma genes  - higher on average in pCR
#   * proliferation (E2F / G2M) genes  - higher on average in No pCR
#   * oestrogen-response genes         - higher on average in No pCR
# Genes are real human Ensembl IDs (the xCell signature genes plus a random
# set of other genes) so that symbol mapping, xCell and enrichment run as they
# would on GEO data. The count matrix columns are labelled with sample titles,
# not accessions, to exercise the matching logic in 01_input_validation.R.
#
# Usage (from the repository root): Rscript example/simulate_example_data.R
# =============================================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(xCell)
})

set.seed(20260115)
out_dir <- "example/data_raw"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Gene universe ---------------------------------------------------------------
immune <- c("CXCL9", "CXCL10", "CXCL11", "CXCL13", "IDO1", "STAT1", "GBP1", "GBP4", "GBP5",
            "IRF1", "CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "CCL4", "CCL5",
            "LAG3", "PDCD1", "CD274", "TAP1", "TAP2", "PSMB8", "PSMB9", "HLA-DRA", "HLA-DPA1",
            "HLA-DPB1", "HLA-DQA1", "CD74", "IFNG", "CD2", "CD3D", "CD3E", "CD3G", "CD247",
            "IL2RG", "TNFRSF9", "ICOS", "CTLA4", "CD27", "ISG15", "IFI44L", "MX1", "OAS1",
            "IFIT1", "IFIT3", "IRF7", "STAT2", "CD38", "MS4A1", "CD79A", "CD19", "IGKC",
            "PTPRC", "LCK", "ZAP70", "CIITA", "SOCS1", "IL10RA", "CD48", "SLAMF7", "CCR5",
            "CXCR3", "FASLG", "KLRD1", "KLRK1", "TIGIT", "CD69", "CD52", "LTB", "SELL",
            "CD86", "CD80", "FCGR3A", "C1QA", "C1QB", "TNF", "IL2RB", "BATF2", "UBE2L6",
            "PARP14", "SAMD9L", "XAF1", "IFI27", "IFI35", "BST2", "HLA-A", "HLA-B", "HLA-C")
prolif <- c("MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "BUB1", "BUB1B", "AURKA", "AURKB",
            "PLK1", "CDC20", "CENPF", "CENPA", "E2F1", "E2F2", "MCM2", "MCM3", "MCM4", "MCM5",
            "MCM6", "MCM7", "PCNA", "TYMS", "RRM2", "TK1", "BIRC5", "KIF2C", "KIF20A", "KIF4A",
            "CCNA2", "CCNE1", "CDKN3", "FOXM1", "MYBL2", "UBE2C", "PTTG1", "CDC6", "CDC25A",
            "CHEK1", "DLGAP5", "HMMR", "NUSAP1", "PRC1", "TPX2", "ZWINT", "MAD2L1", "TTK",
            "NDC80", "NUF2", "ESPL1", "RAD51", "BRCA1", "FEN1", "GINS2", "DHFR", "SMC4",
            "STMN1", "H2AZ1", "HMGB2", "TUBB", "CKS2", "CKS1B", "MELK", "EXO1", "ORC1")
estrogen <- c("ESR1", "PGR", "GATA3", "TFF1", "TFF3", "GREB1", "XBP1", "FOXA1", "AGR2",
              "MYB", "CA12", "SLC39A6", "NAT1", "AR", "KRT18", "KRT8", "SCUBE2", "BCL2",
              "IGF1R", "MLPH", "TBC1D9", "ERBB4", "PDZK1", "RARA", "CCND1", "STC2", "IL6ST")

data("xCell.data", package = "xCell", envir = environment())
sym_all <- unique(c(xCell.data$genes, immune, prolif, estrogen))
ens <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = sym_all, column = "ENSEMBL",
                             keytype = "SYMBOL", multiVals = "first")
genes <- data.frame(symbol = sym_all, gene_id = unname(ens), stringsAsFactors = FALSE)
genes <- genes[!is.na(genes$gene_id) & !duplicated(genes$gene_id), ]

# Background genes not in the xCell signatures.
bg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, "ENSEMBL"),
                            columns = "SYMBOL", keytype = "ENSEMBL")
bg <- bg[!is.na(bg$SYMBOL) & !duplicated(bg$ENSEMBL) & !duplicated(bg$SYMBOL) &
           !bg$SYMBOL %in% genes$symbol & !bg$ENSEMBL %in% genes$gene_id, ]
bg <- bg[sample(nrow(bg), 1500), ]
genes <- rbind(genes, data.frame(symbol = bg$SYMBOL, gene_id = bg$ENSEMBL))
n_genes <- nrow(genes)

# ---- Samples -------------------------------------------------------------------------
n_pcr <- 7; n_no <- 10; n_na <- 1
n <- n_pcr + n_no + n_na
response <- c(rep("pCR", n_pcr), rep("No pCR", n_no), rep(NA, n_na))
ord <- sample(n)
response <- response[ord]
title <- sprintf("Patient_%02d", seq_len(n))
accession <- sprintf("SIM%04d", seq_len(n))
er <- ifelse(response %in% "pCR", sample(c("positive", "negative"), n, TRUE, c(0.3, 0.7)),
             sample(c("positive", "negative"), n, TRUE, c(0.7, 0.3)))
her2 <- sample(c("positive", "negative"), n, TRUE, c(0.25, 0.75))
pr <- ifelse(er == "positive", sample(c("positive", "negative"), n, TRUE, c(0.7, 0.3)), "negative")

# Per-sample module activities (log2 scale).
is_pcr <- response %in% "pCR"
immune_score <- rnorm(n, ifelse(is_pcr, 1.0, -0.3), 0.55)
prolif_score <- rnorm(n, ifelse(is_pcr, -0.5, 0.4), 0.5)
er_score <- ifelse(er == "positive", 1, -1) * 0.8 + rnorm(n, 0, 0.3)

# ---- Gene parameters --------------------------------------------------------------
mu <- exp(rnorm(n_genes, log(150), 1.8))
mu <- pmin(pmax(mu, 0.5), 2e5)
disp <- (0.04 + 1.5 / mu) * exp(rnorm(n_genes, 0, 0.4))
loading <- matrix(0, n_genes, 3)
loading[genes$symbol %in% immune, 1] <- runif(sum(genes$symbol %in% immune), 0.8, 1.6)
loading[genes$symbol %in% prolif, 2] <- runif(sum(genes$symbol %in% prolif), 0.8, 1.4)
loading[genes$symbol %in% estrogen, 3] <- runif(sum(genes$symbol %in% estrogen), 0.6, 1.2)
# Weaker co-regulation of the remaining xCell signature genes, so that cell
# type scores vary between samples.
other_sig <- genes$symbol %in% xCell.data$genes & rowSums(loading) == 0
loading[other_sig, 1] <- rnorm(sum(other_sig), 0, 0.15)

size_factor <- exp(rnorm(n, 0, 0.25))
log2_effect <- loading %*% rbind(immune_score, prolif_score, er_score)
mu_mat <- outer(mu, size_factor) * 2^log2_effect
counts <- matrix(rnbinom(length(mu_mat), mu = mu_mat, size = 1 / disp), n_genes, n)
colnames(counts) <- title

# Mimic GEO formatting: Ensembl version suffix on gene IDs.
gene_id_versioned <- paste0(genes$gene_id, ".", sample(1:20, n_genes, replace = TRUE))

processed <- data.frame(gene_id = gene_id_versioned, gene_symbol = genes$symbol, counts,
                        check.names = FALSE)
processed <- processed[order(processed$gene_id), ]
f_counts <- file.path(out_dir, "example_processed_data.txt.gz")
con <- gzfile(f_counts, "w")
utils::write.table(processed, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

pdata <- data.frame(
  title = title,
  geo_accession = accession,
  source_name_ch1 = "breast tumour biopsy (simulated)",
  `response:ch1` = ifelse(is.na(response), "NA", response),
  `er status:ch1` = er,
  `pr status:ch1` = pr,
  `her2 status:ch1` = her2,
  check.names = FALSE
)
utils::write.table(pdata, file.path(out_dir, "example_pdata.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)

cat("Simulated", n_genes, "genes x", n, "samples\n")
print(table(response, useNA = "ifany"))

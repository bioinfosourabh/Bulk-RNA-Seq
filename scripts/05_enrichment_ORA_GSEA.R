#!/usr/bin/env Rscript
# =============================================================================
# 05_enrichment_ORA_GSEA.R
#
# Pathway enrichment of the DESeq2 results.
#   ORA : hypergeometric test of genes higher / lower in the treatment group
#         (padj < enrichment$padj_cutoff) against the universe of tested genes.
#   GSEA: pre-ranked GSEA (fgsea via clusterProfiler) on all tested genes,
#         ranked by the Wald statistic (default) or shrunken log2 fold change.
# Gene-set collections: MSigDB Hallmark (msigdbr), GO Biological Process
# (org.Hs.eg.db / GO.db) and Reactome (ReactomePA / reactome.db).
#
# Input : results/tables/DESeq2_results_all.tsv
# Output: results/tables/DESeq2_results_annotated.tsv
#         results/tables/enrich_<db>_ORA_UP.tsv, enrich_<db>_ORA_DOWN.tsv,
#         results/tables/enrich_<db>_GSEA.tsv   (<db> = hallmark, GO_BP, reactome)
#         results/tables/enrichment_summary.tsv
#         results/rds/enrichment_results.rds
# =============================================================================

source("scripts/utils.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

cfg <- load_config()
log_step("05. Pathway enrichment (ORA and GSEA)")
ecfg <- cfg$enrichment
dbs <- unlist(ecfg$databases)

require_files(p_table(cfg, "DESeq2_results_all.tsv"))
res <- read_tsv_base(p_table(cfg, "DESeq2_results_all.tsv"))

# ---- Map to Entrez IDs (one Entrez ID per gene) --------------------------------
clean_ids <- strip_ensembl_version(res$gene_id)
if (mean(grepl("^ENSG", clean_ids)) > 0.5) {
  res$ENTREZID <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = clean_ids,
                                               column = "ENTREZID", keytype = "ENSEMBL",
                                               multiVals = "first"))
} else {
  res$ENTREZID <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = res$symbol,
                                               column = "ENTREZID", keytype = "SYMBOL",
                                               multiVals = "first"))
}
log_msg("Genes mapped to Entrez: ", sum(!is.na(res$ENTREZID)), " of ", nrow(res))
write_tsv_base(res, p_table(cfg, "DESeq2_results_annotated.tsv"))

tested <- res[!is.na(res$ENTREZID) & !is.na(res$padj), ]
universe <- unique(tested$ENTREZID)
is_sig <- tested$padj < ecfg$padj_cutoff & abs(tested$log2FoldChange) > ecfg$lfc_cutoff
gene_up <- unique(tested$ENTREZID[is_sig & tested$log2FoldChange > 0])
gene_down <- unique(tested$ENTREZID[is_sig & tested$log2FoldChange < 0])
log_msg("ORA input: ", length(gene_up), " up, ", length(gene_down),
        " down, universe ", length(universe))

# Ranked list for GSEA: highest absolute statistic per Entrez ID, decreasing.
rank_col <- ecfg$gsea_rank_by
ranked <- res[!is.na(res$ENTREZID) & !is.na(res[[rank_col]]), c("ENTREZID", rank_col)]
ranked <- ranked[order(abs(ranked[[rank_col]]), decreasing = TRUE), ]
ranked <- ranked[!duplicated(ranked$ENTREZID), ]
gene_list <- sort(stats::setNames(ranked[[rank_col]], ranked$ENTREZID), decreasing = TRUE)
log_msg("GSEA input: ", length(gene_list), " genes ranked by ", rank_col)

# ---- Helpers -------------------------------------------------------------------
safe_run <- function(label, expr) {
  out <- tryCatch(expr, error = function(e) {
    log_msg("  ", label, " failed: ", conditionMessage(e))
    NULL
  })
  n <- if (is.null(out)) 0 else nrow(as.data.frame(out))
  log_msg(sprintf("  %-22s %d terms", label, n))
  out
}

readable <- function(x) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(x)
  tryCatch(setReadable(x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
           error = function(e) x)
}

write_result <- function(x, name) {
  f <- p_table(cfg, paste0("enrich_", name, ".tsv"))
  df <- if (is.null(x)) data.frame() else as.data.frame(x)
  if (ncol(df) == 0) {
    df <- data.frame(ID = character(), Description = character(), p.adjust = numeric())
  }
  write_tsv_base(df, f)
  nrow(df)
}

ora_args <- list(pAdjustMethod = "BH", pvalueCutoff = ecfg$ora_pvalue_cutoff,
                 qvalueCutoff = ecfg$ora_qvalue_cutoff,
                 minGSSize = ecfg$min_gs_size, maxGSSize = ecfg$max_gs_size)
gsea_args <- list(pAdjustMethod = "BH", pvalueCutoff = ecfg$gsea_pvalue_cutoff,
                  minGSSize = ecfg$min_gs_size, maxGSSize = ecfg$max_gs_size,
                  eps = 0, seed = TRUE, verbose = FALSE)

results <- list()

# ---- MSigDB Hallmark ---------------------------------------------------------
if ("hallmark" %in% dbs) {
  log_msg("Hallmark")
  if (!is.null(ecfg$hallmark_gmt)) {
    t2g <- read.gmt(ecfg$hallmark_gmt)
    colnames(t2g) <- c("term", "gene")
  } else {
    fa <- names(formals(msigdbr::msigdbr))
    h <- if ("collection" %in% fa) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }
    id_col <- intersect(c("ncbi_gene", "entrez_gene"), colnames(h))[1]
    t2g <- data.frame(term = h$gs_name, gene = as.character(h[[id_col]]))
  }
  t2g <- unique(t2g)
  set.seed(cfg$seed)
  results$hallmark_ORA_UP <- safe_run("ORA up", readable(do.call(enricher,
    c(list(gene = gene_up, universe = universe, TERM2GENE = t2g), ora_args))))
  results$hallmark_ORA_DOWN <- safe_run("ORA down", readable(do.call(enricher,
    c(list(gene = gene_down, universe = universe, TERM2GENE = t2g), ora_args))))
  set.seed(cfg$seed)
  results$hallmark_GSEA <- safe_run("GSEA", readable(do.call(GSEA,
    c(list(geneList = gene_list, TERM2GENE = t2g), gsea_args))))
}

# ---- GO Biological Process -----------------------------------------------------
if ("go_bp" %in% dbs) {
  log_msg("GO Biological Process")
  go_args <- list(OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
  results$GO_BP_ORA_UP <- safe_run("ORA up", do.call(enrichGO,
    c(list(gene = gene_up, universe = universe, readable = TRUE), go_args, ora_args)))
  results$GO_BP_ORA_DOWN <- safe_run("ORA down", do.call(enrichGO,
    c(list(gene = gene_down, universe = universe, readable = TRUE), go_args, ora_args)))
  set.seed(cfg$seed)
  results$GO_BP_GSEA <- safe_run("GSEA", readable(do.call(gseGO,
    c(list(geneList = gene_list), go_args, gsea_args))))
}

# ---- Reactome -------------------------------------------------------------------
if ("reactome" %in% dbs) {
  log_msg("Reactome")
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    stop("ReactomePA is required for the Reactome analysis (or remove 'reactome' ",
         "from enrichment$databases).")
  }
  results$reactome_ORA_UP <- safe_run("ORA up", do.call(ReactomePA::enrichPathway,
    c(list(gene = gene_up, universe = universe, organism = "human", readable = TRUE), ora_args)))
  results$reactome_ORA_DOWN <- safe_run("ORA down", do.call(ReactomePA::enrichPathway,
    c(list(gene = gene_down, universe = universe, organism = "human", readable = TRUE), ora_args)))
  set.seed(cfg$seed)
  results$reactome_GSEA <- safe_run("GSEA", readable(do.call(ReactomePA::gsePathway,
    c(list(geneList = gene_list, organism = "human"), gsea_args))))
}

# ---- Export ----------------------------------------------------------------------
summary_df <- data.frame(analysis = names(results),
                         n_terms = vapply(names(results), function(n) write_result(results[[n]], n),
                                          numeric(1)))
print(summary_df, row.names = FALSE)
write_tsv_base(summary_df, p_table(cfg, "enrichment_summary.tsv"))
saveRDS(results, p_rds(cfg, "enrichment_results.rds"))

write_session_info(cfg, "05_enrichment_ORA_GSEA")

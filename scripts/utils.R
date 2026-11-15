# =============================================================================
# utils.R -- shared helpers for the GSE192341 bulk RNA-seq pipeline
#
# Sourced at the top of every numbered script. Provides:
#   * configuration loading (config.yaml or the file given with --config)
#   * output directory creation and path helpers
#   * logging and sessionInfo capture
#   * a PNG + PDF figure writer
#   * Ensembl -> gene symbol mapping
#
# All paths in the configuration are relative to the repository root, and the
# scripts are expected to be run from the repository root, e.g.
#   Rscript scripts/02_deseq2_analysis.R --config config.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

# ---- configuration ----------------------------------------------------------

#' Locate the configuration file.
#' Priority: --config=<file> or --config <file> on the command line, then the
#' PIPELINE_CONFIG environment variable, then ./config.yaml.
get_config_path <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep("^--config", args)
  if (length(hit) > 0) {
    a <- args[hit[1]]
    if (grepl("=", a)) return(sub("^--config=", "", a))
    if (length(args) > hit[1]) return(args[hit[1] + 1])
  }
  env <- Sys.getenv("PIPELINE_CONFIG", unset = "")
  if (nzchar(env)) return(env)
  "config.yaml"
}

#' Read the configuration and derive output paths.
load_config <- function(path = get_config_path()) {
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path,
         "\nRun the scripts from the repository root, or pass --config <file>.")
  }
  cfg <- yaml::read_yaml(path)
  cfg$config_file <- normalizePath(path)

  res <- cfg$paths$results
  cfg$dirs <- list(
    data_raw = cfg$paths$data_raw,
    data     = cfg$paths$data,
    results  = res,
    tables   = file.path(res, "tables"),
    rds      = file.path(res, "rds"),
    png      = file.path(res, "figures", "png"),
    pdf      = file.path(res, "figures", "pdf"),
    logs     = file.path(res, "logs")
  )
  for (d in cfg$dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cfg$seed)) cfg$seed <- 42
  set.seed(cfg$seed)
  cfg
}

# Path helpers
p_data   <- function(cfg, f) file.path(cfg$dirs$data, f)
p_raw    <- function(cfg, f) file.path(cfg$dirs$data_raw, f)
p_table  <- function(cfg, f) file.path(cfg$dirs$tables, f)
p_rds    <- function(cfg, f) file.path(cfg$dirs$rds, f)

#' Stop with a clear message if a required input is missing.
require_files <- function(...) {
  f <- c(...)
  missing <- f[!file.exists(f)]
  if (length(missing) > 0) {
    stop("Required input(s) not found:\n  ", paste(missing, collapse = "\n  "),
         "\nRun the preceding pipeline step(s) first.")
  }
  invisible(TRUE)
}

# ---- logging ----------------------------------------------------------------

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}

log_step <- function(title) {
  cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
}

#' Write sessionInfo() for a script to results/logs/<script>_sessionInfo.txt.
write_session_info <- function(cfg, script_name) {
  out <- file.path(cfg$dirs$logs, paste0(script_name, "_sessionInfo.txt"))
  writeLines(c(
    paste("Script:", script_name),
    paste("Config:", cfg$config_file),
    paste("Seed:", cfg$seed),
    paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    utils::capture.output(utils::sessionInfo())
  ), out)
  log_msg("sessionInfo written to ", out)
}

# ---- table I/O --------------------------------------------------------------

write_tsv_base <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "NA")
}

read_tsv_base <- function(path, ...) {
  utils::read.delim(path, header = TRUE, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "", ...)
}

# ---- figures ------------------------------------------------------------------

#' Save a ggplot (or patchwork) object as PNG and PDF.
save_figure <- function(cfg, plot, name, width, height) {
  force(plot)  # evaluate before a device is opened, so errors leave no empty file
  dpi <- if (is.null(cfg$figures$dpi)) 300 else cfg$figures$dpi
  png_f <- file.path(cfg$dirs$png, paste0(name, ".png"))
  pdf_f <- file.path(cfg$dirs$pdf, paste0(name, ".pdf"))
  ggplot2::ggsave(png_f, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(pdf_f, plot, width = width, height = height, device = grDevices::pdf)
  log_msg("Saved ", name)
  invisible(png_f)
}

# ---- gene identifiers -------------------------------------------------------

#' Strip Ensembl version suffixes (ENSG00000141510.17 -> ENSG00000141510).
strip_ensembl_version <- function(ids) sub("^(ENS[A-Z]*[0-9]+)\\.[0-9]+$", "\\1", ids)

#' Map gene identifiers to HGNC symbols.
#' Uses the GEO-provided annotation (data/gene_annotation.csv) when available,
#' falling back to org.Hs.eg.db; unmapped identifiers are returned unchanged.
map_to_symbol <- function(ids, cfg = NULL) {
  ids <- as.character(ids)
  symbols <- rep(NA_character_, length(ids))

  ann_file <- if (!is.null(cfg)) p_data(cfg, "gene_annotation.csv") else ""
  if (nzchar(ann_file) && file.exists(ann_file)) {
    ann <- utils::read.csv(ann_file, stringsAsFactors = FALSE)
    if (all(c("gene_id", "gene_symbol") %in% colnames(ann))) {
      lookup <- stats::setNames(ann$gene_symbol, ann$gene_id)
      symbols <- unname(lookup[ids])
    }
  }

  todo <- is.na(symbols) | symbols == ""
  clean <- strip_ensembl_version(ids)
  if (any(todo & grepl("^ENSG", clean)) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    keys <- unique(clean[todo & grepl("^ENSG", clean)])
    mapped <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = keys, column = "SYMBOL",
      keytype = "ENSEMBL", multiVals = "first"))
    symbols[todo] <- unname(mapped[clean[todo]])
  }

  still <- is.na(symbols) | symbols == ""
  symbols[still] <- clean[still]
  symbols
}

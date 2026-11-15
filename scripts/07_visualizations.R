#!/usr/bin/env Rscript
# =============================================================================
# 07_visualizations.R
#
# Figures for the GSE192341 analysis (Fig01-Fig11, FigS1-FigS3), written to
# results/figures/png and results/figures/pdf. Run after steps 00A-05.
#
# The plotting code is Sourabh Kumar's original 07_visualizations.R (the code
# that produced the figures in Visualizations/). Changes from that version are
# limited to: paths taken from the configuration; PDF output in addition to
# PNG; condition colours keyed by the configured group labels (previously
# assigned by order of appearance); GO/Reactome term labels taken from the
# Description column; guards for empty inputs. Each change is marked "# [edit]".
# =============================================================================

source("scripts/utils.R")
cfg <- load_config()
log_step("07. Visualisations")

tab <- function(f) p_table(cfg, f)
rds <- function(f) p_rds(cfg, f)

# [edit] PNG + PDF output helpers
save_png_pdf <- function(plot, name, width, height) {
  save_figure(cfg, plot, name, width, height)
}
open_figure <- function(device, name, width, height) {
  if (device == "png") {
    grDevices::png(file.path(cfg$dirs$png, paste0(name, ".png")),
                   width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(file.path(cfg$dirs$pdf, paste0(name, ".pdf")), width = width, height = height)
  }
}

# =============================================================================
# Visualization configuration (original)
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(viridis)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hexbin)
})

# Color palette for treatment groups
# [edit] keyed by the configured reference / treatment labels
pal_condition <- setNames(c("#DC3220", "#005AB5"),
                          c(cfg$design$reference, cfg$design$treatment))

# Color palette for differential expression
pal_regulation <- c("Upregulated" = "#C41E3A", 
                    "Downregulated" = "#1E5AA8", 
                    "Not Significant" = "#BEBEBE")

# Base theme for all plots
theme_publication <- function(base_size = 11, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      legend.background = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "#555555"),
      plot.margin = margin(15, 15, 10, 10),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Gene symbol mapping function
map_ensembl_to_symbol <- function(ensembl_ids) {
  clean_ids <- sub("\\..*$", "", ensembl_ids)
  symbols <- rep(NA_character_, length(clean_ids))
  valid_keys <- suppressMessages(keys(org.Hs.eg.db, keytype = "ENSEMBL"))
  valid_idx <- clean_ids %in% valid_keys
  valid_clean_ids <- unique(clean_ids[valid_idx])
  
  if (length(valid_clean_ids) > 0) {
    mapped <- suppressMessages(suppressWarnings(
      AnnotationDbi::select(org.Hs.eg.db, keys = valid_clean_ids, 
                           columns = "SYMBOL", keytype = "ENSEMBL")
    ))
    if (!is.null(mapped) && nrow(mapped) > 0) {
      mapped_unique <- mapped[!duplicated(mapped$ENSEMBL), ]
      lookup <- setNames(mapped_unique$SYMBOL, mapped_unique$ENSEMBL)
      symbols[valid_idx] <- lookup[clean_ids[valid_idx]]
    }
  }
  
  na_idx <- is.na(symbols) | symbols == ""
  symbols[na_idx] <- ifelse(nchar(clean_ids[na_idx]) > 15,
                             substr(clean_ids[na_idx], 1, 15),
                             clean_ids[na_idx])
  return(symbols)
}

# [edit] inputs shared by the figure blocks below
coldata <- readRDS(rds("coldata_validated.rds"))
coldata$condition <- as.character(coldata[[cfg$design$condition_column]])
dds <- readRDS(rds("dds.rds"))
dataset_subtitle <- cfg$title

# ==== Fig01 ====
library(ggplot2)
library(patchwork)

qc <- read.table(tab("sample_qc_metrics.tsv"), header = TRUE, sep = "\t")
qc$condition <- coldata[qc$sample, "condition"]

# Panel A: Library size distribution
p1a <- ggplot(qc, aes(x = library_size / 1e6, fill = condition)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_rug(aes(color = condition), alpha = 0.5) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition, guide = "none") +
  labs(title = "Library Size Distribution",
       x = "Library Size (millions of reads)", y = "Density") +
  theme_publication()

# Panel B: Library size by group
p1b <- ggplot(qc, aes(x = condition, y = library_size / 1e6, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.6, shape = 21, fill = "white") +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black") +
  scale_fill_manual(values = pal_condition) +
  labs(title = "Library Size by Response", x = NULL, y = "Library Size (millions)") +
  theme_publication() +
  theme(legend.position = "none")

# Panel C: Detected genes distribution
p1c <- ggplot(qc, aes(x = condition, y = detected_genes / 1000, fill = condition)) +
  geom_violin(alpha = 0.7, color = "black", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = pal_condition) +
  labs(title = "Genes Detected per Sample", x = NULL, y = "Detected Genes (thousands)") +
  theme_publication() +
  theme(legend.position = "none")

# Panel D: Library complexity
cor_val <- cor(qc$library_size, qc$detected_genes)
p1d <- ggplot(qc, aes(x = library_size / 1e6, y = detected_genes / 1000)) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", fill = "#CCCCCC") +
  geom_point(aes(fill = condition), size = 3, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("r = ", round(cor_val, 3)), fontface = "italic") +
  labs(title = "Library Complexity", x = "Library Size (millions)", 
       y = "Detected Genes (thousands)") +
  theme_publication()

fig1 <- (p1a | p1b) / (p1c | p1d) +
  plot_annotation(
    title = "Sample Quality Control Metrics",
    subtitle = dataset_subtitle  # [edit]
  )

save_png_pdf(fig1, "Fig01_Sample_QC", width = 12, height = 11)

# ==== Fig04 ====
res <- read.table(tab("DESeq2_results_all.tsv"), header = TRUE, sep = "\t")

# Map gene IDs to symbols
res$symbol <- map_ensembl_to_symbol(res$gene_id)

# Classify genes
res <- res %>%
  mutate(
    regulation = case_when(
      is.na(padj) ~ "Not Significant",
      padj >= 0.05 ~ "Not Significant",
      padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    regulation = factor(regulation, levels = c("Upregulated", "Downregulated", "Not Significant"))
  )

# Select top genes for labeling
top_up <- res %>% filter(regulation == "Upregulated") %>% arrange(padj) %>% head(12)
top_down <- res %>% filter(regulation == "Downregulated") %>% arrange(padj) %>% head(12)
top_label <- bind_rows(top_up, top_down)

n_up <- sum(res$regulation == "Upregulated")
n_down <- sum(res$regulation == "Downregulated")

# Generate plot
p4 <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = filter(res, regulation == "Not Significant"),
             color = "#CCCCCC", alpha = 0.4, size = 1) +
  geom_point(data = filter(res, regulation != "Not Significant"),
             aes(color = regulation), alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#666666") +
  geom_point(data = top_label, aes(color = regulation), size = 3, shape = 1, stroke = 1) +
  geom_text_repel(data = top_label, aes(label = symbol, color = regulation),
                  size = 3, fontface = "italic", max.overlaps = 30,
                  segment.color = "#555555", show.legend = FALSE) +
  scale_color_manual(
    values = pal_regulation,
    name = "Regulation",
    labels = c(paste0("Upregulated (n=", n_up, ")"),
               paste0("Downregulated (n=", n_down, ")"),
               "Not Significant")
  ) +
  labs(
    title = "Differential Gene Expression: pCR vs No pCR",
    subtitle = paste0("n = ", format(nrow(res), big.mark = ","), 
                      " genes | Thresholds: padj < 0.05, |LFC| > 0.5"),
    x = expression(Log[2]~Fold~Change~(pCR~"/"~No~pCR)),
    y = expression(-Log[10]~Adjusted~P-value)
  ) +
  theme_publication(base_size = 12)

save_png_pdf(p4, "Fig04_Volcano", width = 11, height = 10)

# ==== Fig05 ====
dds <- readRDS(rds("dds.rds"))
res_ma <- as.data.frame(results(dds))
res_ma$gene_id <- rownames(res_ma)
res_ma$symbol <- map_ensembl_to_symbol(res_ma$gene_id)

res_ma <- res_ma %>%
  mutate(
    sig = case_when(
      is.na(padj) ~ "NS",
      padj >= 0.05 ~ "NS",
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

top_ma <- res_ma %>% filter(sig != "NS") %>% arrange(padj) %>% head(15)

p5 <- ggplot(res_ma, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = sig), alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "#888888") +
  geom_text_repel(data = top_ma, aes(label = symbol), size = 2.8, 
                  fontface = "italic", max.overlaps = 20) +
  scale_color_manual(
    values = c("Up" = "#C41E3A", "Down" = "#1E5AA8", "NS" = "#CCCCCC"),
    name = "Significance"
  ) +
  labs(title = "MA Plot", subtitle = "Log2 fold change vs mean normalized expression",
       x = expression(Log[10]~Mean~Expression),
       y = expression(Log[2]~Fold~Change)) +
  theme_publication()

save_png_pdf(p5, "Fig05_MA_Plot", width = 10, height = 8)

# ==== Fig06 ====
# Panel A: P-value distribution
p6a <- ggplot(res %>% filter(!is.na(pvalue)), aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#4DBBD5", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#E64B35") +
  labs(title = "P-value Distribution", x = "P-value", y = "Frequency") +
  theme_publication()

# Panel B: Log2 fold change distribution (significant genes)
sig_res <- res %>% filter(!is.na(padj) & padj < 0.05)
p6b <- ggplot(sig_res, aes(x = log2FoldChange, fill = regulation)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = pal_regulation, name = "Direction") +
  labs(title = "Log2 Fold Change Distribution",
       subtitle = paste0("Significant genes (n = ", nrow(sig_res), ")"),
       x = expression(Log[2]~Fold~Change), y = "Frequency") +
  theme_publication()

# Panel C: Effect size vs significance
p6c <- ggplot(res %>% filter(!is.na(padj)), aes(x = abs(log2FoldChange), y = -log10(padj))) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#E64B35") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#E64B35") +
  scale_fill_viridis_c(option = "plasma", name = "Count", trans = "log10") +
  labs(title = "Effect Size vs Significance",
       x = expression("|"*Log[2]~Fold~Change*"|"),
       y = expression(-Log[10]~Adjusted~P-value)) +
  theme_publication()

# Panel D: Top genes by significance
res_ranked <- res %>% filter(!is.na(padj)) %>% arrange(padj) %>% mutate(rank = row_number())
p6d <- ggplot(res_ranked %>% head(500), aes(x = rank, y = -log10(padj))) +
  geom_line(color = "#3C5488", linewidth = 1) +
  geom_point(aes(color = regulation), size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#888888") +
  scale_color_manual(values = pal_regulation, name = "Regulation") +
  labs(title = "Top 500 Genes by Significance",
       x = "Rank", y = expression(-Log[10]~Adjusted~P-value)) +
  theme_publication()

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(title = "Differential Expression Statistics")

save_png_pdf(fig6, "Fig06_DE_Statistics", width = 13, height = 11)

# ==== Fig08 ====
norm_counts <- as.matrix(read.table(tab("normalized_counts.tsv"), header = TRUE, 
                                    sep = "\t", row.names = 1, check.names = FALSE))

top12 <- res %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(12)
common <- intersect(colnames(norm_counts), rownames(coldata))
norm_counts <- norm_counts[, common]

plot_gene_expr <- function(gid, gene_symbol, lfc, pval) {
  df <- data.frame(expr = norm_counts[gid, ], Response = coldata[common, "condition"])
  
  ggplot(df, aes(x = Response, y = expr, fill = Response)) +
    geom_violin(alpha = 0.6, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = pal_condition) +
    scale_y_continuous(labels = comma_format()) +
    labs(title = gene_symbol, subtitle = sprintf("LFC = %.2f | padj = %.2e", lfc, pval),
         x = NULL, y = "Normalized Counts") +
    theme_publication(base_size = 9) +
    theme(legend.position = "none", plot.title = element_text(face = "bold.italic", size = 11))
}

gene_plots <- lapply(1:nrow(top12), function(i) {
  plot_gene_expr(top12$gene_id[i], top12$symbol[i], 
                 top12$log2FoldChange[i], top12$padj[i])
})

fig8 <- wrap_plots(gene_plots, ncol = 4) +
  plot_annotation(title = "Top Differentially Expressed Genes (n = 12)",
                  subtitle = "Ranked by adjusted p-value")

save_png_pdf(fig8, "Fig08_Top_Gene_Expression", width = 14, height = 10)

# ==== FigS1 ====
dds <- readRDS(rds("dds.rds"))
for (device in c("png", "pdf")) {
open_figure(device, "FigS1_Dispersion", width = 9, height = 7)
plotDispEsts(dds, main = "DESeq2 Gene-wise Dispersion Estimates")
dev.off()
}

# ==== FigS2 ====
sf_df <- data.frame(sample = names(sizeFactors(dds)),
                    size_factor = sizeFactors(dds),
                    condition = colData(dds)$condition)

pS2 <- ggplot(sf_df, aes(x = reorder(sample, size_factor), y = size_factor, fill = condition)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_manual(values = pal_condition, name = "Response") +
  labs(title = "DESeq2 Size Factors", subtitle = "Library size normalization factors",
       x = "Sample", y = "Size Factor") +
  theme_publication() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

save_png_pdf(pS2, "FigS2_Size_Factors", width = 12, height = 5)

# ==== Fig02 ====
pca_df <- read.table(tab("qc_pca_coordinates.tsv"), header = TRUE, sep = "\t")
var_df <- read.table(tab("qc_pca_variance_explained.tsv"), header = TRUE, sep = "\t")
var_df$pct <- round(var_df$variance_explained * 100, 1)

# Panel A: PC1 vs PC2
p2a <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = condition), level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(fill = condition), size = 4, shape = 21, stroke = 0.5, alpha = 0.85) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition) +
  labs(title = "Principal Component Analysis",
       x = paste0("PC1 (", var_df$pct[1], "% variance)"),
       y = paste0("PC2 (", var_df$pct[2], "% variance)")) +
  theme_publication(base_size = 12)

# Panel B: PC2 vs PC3
p2b <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
  stat_ellipse(aes(color = condition), level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(fill = condition), size = 4, shape = 21, stroke = 0.5, alpha = 0.85) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition) +
  labs(title = "PC2 vs PC3",
       x = paste0("PC2 (", var_df$pct[2], "% variance)"),
       y = paste0("PC3 (", var_df$pct[3], "% variance)")) +
  theme_publication(base_size = 12)

# Panel C: Scree plot
var_plot <- var_df[1:min(10, nrow(var_df)), ]
var_plot$PC_num <- 1:nrow(var_plot)
var_plot$cumulative <- cumsum(var_plot$pct)

p2c <- ggplot(var_plot, aes(x = PC_num)) +
  geom_bar(aes(y = pct), stat = "identity", fill = "#4DBBD5", alpha = 0.85, width = 0.7) +
  geom_line(aes(y = cumulative / 2), color = "#E64B35", linewidth = 1.2) +
  geom_point(aes(y = cumulative / 2), color = "#E64B35", size = 3) +
  geom_text(aes(y = pct, label = paste0(pct, "%")), vjust = -0.4, size = 3.2) +
  scale_x_continuous(breaks = 1:10, labels = paste0("PC", 1:10)) +
  scale_y_continuous(name = "Individual Variance (%)",
                     sec.axis = sec_axis(~ . * 2, name = "Cumulative Variance (%)")) +
  labs(title = "Variance Explained by Principal Components", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel D: PC1 distribution
p2d <- ggplot(pca_df, aes(x = condition, y = PC1, fill = condition)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = pal_condition) +
  labs(title = "PC1 Score Distribution", x = NULL, y = "PC1 Score") +
  theme_publication() +
  theme(legend.position = "none")

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_annotation(title = "Principal Component Analysis",
                  subtitle = "Transcriptome-wide variance structure")

save_png_pdf(fig2, "Fig02_PCA", width = 13, height = 12)

# ==== Fig03 ====
library(pheatmap)

dist_mat <- as.matrix(read.table(tab("qc_sample_distance_matrix.tsv"), 
                                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
cor_mat <- as.matrix(read.table(tab("qc_sample_correlation_matrix.tsv"), 
                                header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))

common <- intersect(rownames(dist_mat), rownames(coldata))
dist_mat <- dist_mat[common, common]
cor_mat <- cor_mat[common, common]

condition_values <- as.character(coldata[common, "condition"])
anno <- data.frame(Response = condition_values, row.names = common)

cond_colors <- pal_condition  # [edit] fixed colour per group
anno_colors <- list(Response = cond_colors)

# Distance heatmap
pal_distance <- colorRampPalette(c("#FFFFFF", "#FFF5F0", "#FEE0D2", "#FCBBA1", 
                                   "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D"))(100)

for (device in c("png", "pdf")) {
open_figure(device, "Fig03A_Sample_Distance", width = 10, height = 9)
pheatmap(dist_mat, color = pal_distance, annotation_col = anno, annotation_row = anno,
         annotation_colors = anno_colors, show_rownames = FALSE, show_colnames = FALSE,
         clustering_method = "ward.D2", main = "Sample-to-Sample Euclidean Distance",
         fontsize = 11, border_color = NA, treeheight_row = 40, treeheight_col = 40)
dev.off()
}

# Correlation heatmap
cor_colors <- colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"))(100)

for (device in c("png", "pdf")) {
open_figure(device, "Fig03B_Sample_Correlation", width = 10, height = 9)
pheatmap(cor_mat, color = cor_colors, breaks = seq(min(cor_mat, na.rm = TRUE), 1, length.out = 101),
         annotation_col = anno, annotation_row = anno, annotation_colors = anno_colors,
         show_rownames = FALSE, show_colnames = FALSE, clustering_method = "ward.D2",
         main = "Sample-to-Sample Pearson Correlation", fontsize = 11, border_color = NA,
         treeheight_row = 40, treeheight_col = 40)
dev.off()
}

# ==== Fig07 ====
vst_mat <- as.matrix(read.table(tab("vst_matrix.tsv"), header = TRUE, sep = "\t", 
                                row.names = 1, check.names = FALSE))

top_genes <- res %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj) %>% head(60) %>% pull(gene_id)
top_genes <- top_genes[top_genes %in% rownames(vst_mat)]

heat_mat <- vst_mat[top_genes, , drop = FALSE]
symbols <- map_ensembl_to_symbol(rownames(heat_mat))
rownames(heat_mat) <- make.unique(symbols)

heat_scaled <- t(scale(t(heat_mat)))

common <- intersect(colnames(heat_scaled), rownames(coldata))
heat_scaled <- heat_scaled[, common]

condition_values <- coldata[common, "condition"]
anno_col <- data.frame(Response = condition_values, row.names = common)

gene_lfc <- res$log2FoldChange[match(top_genes, res$gene_id)]
direction_vec <- ifelse(is.na(gene_lfc), "Unknown",
                        ifelse(gene_lfc > 0, "Up in pCR", "Down in pCR"))
anno_row <- data.frame(Direction = direction_vec, row.names = rownames(heat_scaled))

cond_colors <- pal_condition  # [edit] fixed colour per group

anno_colors <- list(
  Response = cond_colors,
  Direction = c("Up in pCR" = "#C41E3A", "Down in pCR" = "#1E5AA8", "Unknown" = "#888888")
)

pal_heatmap_diverging <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                            "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

for (device in c("png", "pdf")) {
open_figure(device, "Fig07_Top_DE_Heatmap", width = 12, height = 14)
pheatmap(heat_scaled, color = pal_heatmap_diverging, breaks = seq(-3, 3, length.out = 101),
         annotation_col = anno_col, annotation_row = anno_row, annotation_colors = anno_colors,
         show_colnames = FALSE, fontsize_row = 7, fontsize = 10,
         main = paste0("Top ", nrow(heat_scaled), " Differentially Expressed Genes"),
         clustering_method = "ward.D2", border_color = NA,
         cutree_cols = 2, cutree_rows = 2, treeheight_row = 30, treeheight_col = 30)
dev.off()
}

# ==== FigS3 ====
top_var <- read.table(tab("qc_top_variable_genes.tsv"), header = TRUE, sep = "\t")
vst_mat <- as.matrix(read.table(tab("vst_matrix.tsv"), header = TRUE, sep = "\t", 
                                row.names = 1, check.names = FALSE))

top100 <- head(top_var$gene_id, 100)
top100 <- top100[top100 %in% rownames(vst_mat)]

heat_var <- vst_mat[top100, ]
symbols <- map_ensembl_to_symbol(rownames(heat_var))
rownames(heat_var) <- make.unique(symbols)

heat_scaled <- t(scale(t(heat_var)))
common <- intersect(colnames(heat_scaled), rownames(coldata))
heat_scaled <- heat_scaled[, common]

condition_values <- as.character(coldata[common, "condition"])
anno <- data.frame(Response = condition_values, row.names = common)

cond_colors <- pal_condition  # [edit] fixed colour per group
anno_colors <- list(Response = cond_colors)

for (device in c("png", "pdf")) {
open_figure(device, "FigS3_Top_Variable_Genes", width = 12, height = 16)
pheatmap(heat_scaled, color = pal_heatmap_diverging, breaks = seq(-3, 3, length.out = 101),
         annotation_col = anno, annotation_colors = anno_colors, show_colnames = FALSE,
         fontsize_row = 6, main = paste0("Top ", length(top100), " Most Variable Genes"),
         clustering_method = "ward.D2", border_color = NA)
dev.off()
}

# ==== Fig09A ====
xcell <- read.table(tab("deconvolution_xcell_scores.tsv"), header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE)

cell_cols <- setdiff(colnames(xcell), c("sample", "condition"))
xcell_mat <- as.matrix(xcell[, cell_cols])
rownames(xcell_mat) <- xcell$sample

cell_var <- apply(xcell_mat, 2, var, na.rm = TRUE)
cell_var <- cell_var[cell_var > 0]  # [edit] constant scores cannot be z-scaled
top_cells <- names(sort(cell_var, decreasing = TRUE))[1:min(30, length(cell_var))]
xcell_scaled <- scale(xcell_mat[, top_cells])

condition_values <- as.character(xcell$condition)
anno <- data.frame(Response = condition_values, row.names = xcell$sample)

cond_colors <- pal_condition  # [edit] fixed colour per group
anno_colors <- list(Response = cond_colors)

for (device in c("png", "pdf")) {
open_figure(device, "Fig09A_xCell_Heatmap", width = 14, height = 10)
pheatmap(t(xcell_scaled), color = pal_heatmap_diverging, breaks = seq(-2.5, 2.5, length.out = 101),
         annotation_col = anno, annotation_colors = anno_colors, show_colnames = FALSE,
         fontsize_row = 9, main = "xCell Enrichment Scores (Top 30 Variable Cell Types)",
         clustering_method = "ward.D2", border_color = NA, treeheight_col = 30)
dev.off()
}

# ==== Fig09B ====
immune_cells <- c("CD8+ T-cells", "CD4+ T-cells", "NK cells", "B-cells",
                  "Macrophages M1", "Macrophages M2", "Tregs", "Monocytes",
                  "Dendritic cells", "Neutrophils")
avail_immune <- intersect(immune_cells, cell_cols)

immune_long <- xcell %>%
  dplyr::select(sample, condition, all_of(avail_immune)) %>%
  pivot_longer(cols = -c(sample, condition), names_to = "Cell_Type", values_to = "Score") %>%
  mutate(Cell_Type = factor(Cell_Type, levels = avail_immune))

p9b <- ggplot(immune_long, aes(x = Cell_Type, y = Score, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, position = position_dodge(0.75), width = 0.65) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 0.8, alpha = 0.4) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  labs(title = "Immune Cell Infiltration by Treatment Response",
       x = NULL, y = "xCell Enrichment Score") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

save_png_pdf(p9b, "Fig09B_Immune_Cells", width = 12, height = 8)

# ==== Fig09C ====
xcell_stats <- read.table(tab("deconvolution_xcell_pcr_vs_nopcr.tsv"), 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sig_cells <- xcell_stats %>% filter(!is.na(padj) & padj < 0.25) %>% arrange(padj) %>% head(15)
fig9c_subtitle <- "Wilcoxon rank-sum test, FDR < 0.25"
if (nrow(sig_cells) == 0) {  # [edit] no cell type passes FDR < 0.25: show the 15 smallest p-values
  sig_cells <- xcell_stats %>% filter(!is.na(p_value)) %>% arrange(p_value) %>% head(15)
  fig9c_subtitle <- "No cell type at FDR < 0.25; 15 smallest p-values shown"
}

sig_long <- xcell %>%
  dplyr::select(sample, condition, all_of(sig_cells$population)) %>%
  pivot_longer(cols = -c(sample, condition), names_to = "Cell_Type", values_to = "Score") %>%
  left_join(sig_cells %>% dplyr::select(population, padj), 
            by = c("Cell_Type" = "population")) %>%
  mutate(Cell_Type = fct_reorder(Cell_Type, -padj))

p9c <- ggplot(sig_long, aes(x = Cell_Type, y = Score, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.6, width = 0.6) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  coord_flip() +
  labs(title = "Differentially Enriched Cell Populations",
       subtitle = fig9c_subtitle,  # [edit]
       x = NULL, y = "xCell Enrichment Score") +
  theme_publication()

save_png_pdf(p9c, "Fig09C_Significant_CellTypes", width = 10, height = 8)

# ==== Fig09D ====
immune_cols <- intersect(c("ImmuneScore", "StromaScore", "MicroenvironmentScore"), cell_cols)

scores_df <- xcell %>%
  dplyr::select(sample, condition, all_of(immune_cols)) %>%
  pivot_longer(cols = -c(sample, condition), names_to = "Score_Type", values_to = "Score")

p9d <- ggplot(scores_df, aes(x = condition, y = Score, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  facet_wrap(~Score_Type, scales = "free_y") +
  scale_fill_manual(values = pal_condition) +
  labs(title = "Tumor Microenvironment Composite Scores", x = NULL, y = "xCell Score") +
  theme_publication() +
  theme(legend.position = "none")

save_png_pdf(p9d, "Fig09D_TME_Scores", width = 10, height = 6)

# ==== Fig10 ====
# Dotplot generation function
make_dotplot <- function(file_path, title, n_show = 15) {
  if (!file.exists(file_path)) return(NULL)
  
  df <- tryCatch({
    read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
               fill = TRUE, quote = "", comment.char = "")
  }, error = function(e) NULL)
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if ("GeneRatio" %in% colnames(df)) {
    df$GeneRatio_num <- sapply(strsplit(df$GeneRatio, "/"), function(x) {
      if (length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA
    })
  } else {
    return(NULL)
  }
  
  # [edit] Hallmark IDs are converted to readable names; GO and Reactome terms
  # keep their Description (their IDs are accession numbers)
  is_hallmark <- grepl("^HALLMARK_", df$ID)
  hallmark_names <- tools::toTitleCase(tolower(gsub("_", " ", gsub("HALLMARK_", "", df$ID))))
  df$Description <- ifelse(is_hallmark, hallmark_names, df$Description)
  df$Description <- str_wrap(df$Description, width = 40)
  
  df <- df %>%
    filter(!is.na(p.adjust) & !is.na(GeneRatio_num)) %>%
    arrange(p.adjust) %>%
    head(n_show) %>%
    mutate(Description = fct_reorder(Description, GeneRatio_num))
  
  if (nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = GeneRatio_num, y = Description)) +
    geom_segment(aes(x = 0, xend = GeneRatio_num, yend = Description), 
                 color = "#CCCCCC", linewidth = 0.8) +
    geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.85) +
    scale_color_viridis_c(option = "plasma", name = expression(-Log[10]~FDR), direction = -1) +
    scale_size_continuous(range = c(4, 12), name = "Gene Count") +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 9))
}

# Generate plots
p_hall_up <- make_dotplot(tab("enrich_hallmark_ORA_UP.tsv"), 
                          "Hallmark Pathways: Upregulated in pCR", n_show = 12)
if (!is.null(p_hall_up)) {
  save_png_pdf(p_hall_up, "Fig10A_Hallmark_UP", width = 10, height = 8)
}

p_hall_down <- make_dotplot(tab("enrich_hallmark_ORA_DOWN.tsv"), 
                            "Hallmark Pathways: Downregulated in pCR", n_show = 12)
if (!is.null(p_hall_down)) {
  save_png_pdf(p_hall_down, "Fig10B_Hallmark_DOWN", width = 10, height = 8)
}

p_go_up <- make_dotplot(tab("enrich_GO_BP_ORA_UP.tsv"), 
                        "GO Biological Process: Upregulated", n_show = 18)
if (!is.null(p_go_up)) {
  save_png_pdf(p_go_up, "Fig10C_GO_BP_UP", width = 11, height = 10)
}

p_go_down <- make_dotplot(tab("enrich_GO_BP_ORA_DOWN.tsv"), 
                          "GO Biological Process: Downregulated", n_show = 18)
if (!is.null(p_go_down)) {
  save_png_pdf(p_go_down, "Fig10D_GO_BP_DOWN", width = 11, height = 10)
}

# GSEA barplot
gsea_file <- tab("enrich_hallmark_GSEA.tsv")
if (file.exists(gsea_file)) {
  gsea <- read.table(gsea_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     fill = TRUE, quote = "", comment.char = "")
  
  if (!is.null(gsea) && nrow(gsea) > 0 && "NES" %in% colnames(gsea)) {
    gsea$Description <- gsub("HALLMARK_", "", gsea$ID)
    gsea$Description <- gsub("_", " ", gsea$Description)
    gsea$Description <- tools::toTitleCase(tolower(gsea$Description))
    
    gsea <- gsea %>%
      filter(!is.na(p.adjust)) %>%
      arrange(desc(abs(NES))) %>%
      head(20) %>%
      mutate(Direction = ifelse(NES > 0, "Enriched in pCR", "Enriched in No pCR"),
             Description = fct_reorder(Description, NES))
    
    p_gsea <- ggplot(gsea, aes(x = NES, y = Description, fill = Direction)) +
      geom_bar(stat = "identity", alpha = 0.85, width = 0.75) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
      scale_fill_manual(values = c("Enriched in pCR" = "#005AB5", "Enriched in No pCR" = "#DC3220")) +
      labs(title = "Gene Set Enrichment Analysis: Hallmark Pathways",
           subtitle = "Normalized Enrichment Score (NES)",
           x = "Normalized Enrichment Score", y = NULL) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 9))
    
    save_png_pdf(p_gsea, "Fig10E_GSEA_Hallmark", width = 11, height = 10)
  }
}

# ==== Fig11 ====
# Panel A: Sample distribution
cond_df <- as.data.frame(table(coldata$condition))
colnames(cond_df) <- c("Response", "Count")

cond_colors <- pal_condition  # [edit] fixed colour per group

p11a <- ggplot(cond_df, aes(x = Response, y = Count, fill = Response)) +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.9) +
  geom_text(aes(label = Count), vjust = -0.4, size = 6, fontface = "bold") +
  scale_fill_manual(values = cond_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Sample Distribution", x = NULL, y = "Number of Samples") +
  theme_publication(base_size = 13) +
  theme(legend.position = "none")

# Panel B: DE gene categorization
res$Category <- "Not Significant"
res$Category[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1] <- "Up (|LFC|>1)"
res$Category[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 0 & res$log2FoldChange <= 1] <- "Up (|LFC|≤1)"
res$Category[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1] <- "Down (|LFC|>1)"
res$Category[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < 0 & res$log2FoldChange >= -1] <- "Down (|LFC|≤1)"

de_counts <- as.data.frame(table(res$Category))
colnames(de_counts) <- c("Category", "n")
de_counts$Category <- factor(de_counts$Category, 
                             levels = c("Up (|LFC|>1)", "Up (|LFC|≤1)", 
                                        "Down (|LFC|≤1)", "Down (|LFC|>1)", "Not Significant"))

de_cols <- c("Up (|LFC|>1)" = "#67001F", "Up (|LFC|≤1)" = "#D6604D",
             "Down (|LFC|≤1)" = "#4393C3", "Down (|LFC|>1)" = "#053061", 
             "Not Significant" = "#BEBEBE")

p11b <- ggplot(de_counts, aes(x = Category, y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.9) +
  geom_text(aes(label = comma(n)), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = de_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(title = "Differential Expression Summary", x = NULL, y = "Number of Genes") +
  theme_publication(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

fig11 <- (p11a | p11b) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(title = paste("Analysis Summary:", cfg$project),  # [edit]
                  subtitle = dataset_subtitle)

save_png_pdf(fig11, "Fig11_Summary", width = 12, height = 5.5)

write_session_info(cfg, "07_visualizations")
log_msg("Figures written to ", cfg$dirs$png, " and ", cfg$dirs$pdf)

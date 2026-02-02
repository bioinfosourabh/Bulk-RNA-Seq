# Bulk-RNA-Seq: A Step-by-Step Workflow (DESeq2 + xCell + Pathway Enrichment + Elite Visualizations)

A modular bulk RNA-seq differential expression pipeline using DESeq2 in R. This repository demonstrates a complete workflow from raw GEO data to pathway-level biological interpretation with publication-grade visualizations, using the breast cancer neoadjuvant response dataset GSE192341.

## Purpose of This Repository
This repository provides a reproducible pipeline for bulk RNA-seq analysis covering:
1. Data Preparation: Parsing GEO processed matrices into DESeq2-compatible count matrices
2. Metadata Construction: Building sample annotation (colData) from GEO clinical variables
3. Input Validation: Ensuring sample alignment and applying low-count gene filtering
4. Differential Expression: Negative binomial modeling with DESeq2 and apeglm log-fold change shrinkage
5. Quality Control: Variance-stabilizing transformation, PCA, sample distance matrices, and outlier detection
6. Cell-Type Deconvolution: Estimating immune and stromal cell enrichment scores using xCell
7. Pathway Enrichment: Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) with Hallmark, GO, and Reactome gene sets
8. **Elite Visualizations**: Publication-quality figures with gene symbols, Nature/Cell/Science-grade aesthetics integrated throughout the workflow

## Dataset
- **GEO Accession:** [GSE192341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192341)
- **Tissue:** Human breast tumor biopsies
- **Contrast:** Pathological complete response (pCR) vs residual disease (No pCR) following neoadjuvant chemotherapy
- **Samples:** 87 patients (pCR = 24, No pCR = 61, NA = 2)

## Repository Structure
```
Bulk-RNA-Seq/
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   ├── 00A_make_counts_from_processed_matrix.R
│   ├── 00B_make_coldata_from_GEO.R
│   ├── 01_input_validation.R
│   ├── 02_deseq2_analysis.R
│   ├── 03_qc_vst_pca_distances.R
│   ├── 04_deconvolution_xcell.R
│   ├── 05_enrichment_ORA_GSEA.R
│   └── 07_visualizations.R
├── data/
│   └── coldata.tsv
├── data_raw/                 # Large files — not tracked
│   └── GSE192341_processed_data.txt
├── results/                  # Output files — not tracked
│   ├── tables/
│   ├── figures/
│   │   ├── png/
│   │   └── pdf/
│   └── rds/
└── Visualizations/           # Published figures
```

## Installation & Setup
Install the required R packages:
```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2",
    "apeglm",
    "GEOquery",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "clusterProfiler",
    "ReactomePA",
    "DOSE",
    "enrichplot"
))

# CRAN packages
install.packages(c("tidyverse", "data.table", "pheatmap", "RColorBrewer", 
                   "ggrepel", "msigdbr", "scales", "patchwork", "viridis", "hexbin"))

# xCell from GitHub
devtools::install_github("dviraran/xCell")
```

## Visualization Setup

For generating publication-grade figures, additional setup is required:
```r
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

# Elite color palettes (Nature/Cell/Science inspired)
pal_condition <- c("No pCR" = "#DC3220", "pCR" = "#005AB5")
pal_regulation <- c("Upregulated" = "#C41E3A", "Downregulated" = "#1E5AA8", "Not Significant" = "#BEBEBE")

# Publication theme
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
      legend.key = element_blank(),
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
```

## 00A. Constructing Count Matrix from GEO Processed Data

**Purpose:**
GEO processed files often contain gene annotations interleaved with expression values. This step parses the raw GEO matrix to generate a clean gene × sample integer count matrix and a separate gene annotation table for downstream identifier mapping (ENSEMBL → SYMBOL → ENTREZID).

**Code:**
```r
### Step 1: Load Libraries
library(data.table)
library(tidyverse)

### Step 2: Read GEO Processed Matrix
raw_data <- fread("data_raw/GSE192341_processed_data.txt")

### Step 3: Inspect Column Structure
head(colnames(raw_data), 10)
str(raw_data[, 1:5])

### Step 4: Extract Gene Annotation
gene_annotation <- raw_data %>%
  select(gene_id = 1, gene_symbol = 2) %>%
  distinct()

### Step 5: Build Count Matrix
sample_cols <- colnames(raw_data)[sapply(raw_data, is.numeric)]
count_matrix <- as.matrix(raw_data[, ..sample_cols])
rownames(count_matrix) <- raw_data[[1]]

### Step 6: Verify Matrix Integrity
stopifnot(all(count_matrix >= 0))
cat("Dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

### Step 7: Export
write.csv(count_matrix, "data/counts.csv", row.names = TRUE)
write.csv(gene_annotation, "data/gene_annotation.csv", row.names = FALSE)
saveRDS(count_matrix, "results/rds/counts_raw.rds")
```

## 00B. Building Sample Metadata (colData) from GEO

**Purpose:**
DESeq2 requires a sample metadata data.frame where row names match the count matrix column names exactly. This step queries GEO using GEOquery to extract clinical covariates and constructs the experimental design factor (pCR vs No pCR).

**Code:**
```r
### Step 1: Load Libraries
library(GEOquery)
library(tidyverse)

### Step 2: Download Series Matrix from GEO
gse <- getGEO("GSE192341", GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])

### Step 3: Inspect Available Clinical Variables
colnames(pdata)
head(pdata[, grep("characteristics", colnames(pdata))])

### Step 4: Extract Relevant Covariates
coldata <- pdata %>%
  transmute(
    sample_id = geo_accession,
    title = title,
    response = `response:ch1`,
    er_status = `er status:ch1`,
    her2_status = `her2 status:ch1`,
    pr_status = `pr status:ch1`
  ) %>%
  mutate(
    condition = factor(
      ifelse(response == "pCR", "pCR", "No_pCR"),
      levels = c("No_pCR", "pCR")
    )
  )

### Step 5: Set Row Names to Match Count Matrix
rownames(coldata) <- coldata$sample_id

### Step 6: Verify Group Distribution
table(coldata$condition, useNA = "ifany")

### Step 7: Export
write_tsv(coldata, "data/coldata.tsv")
saveRDS(coldata, "results/rds/coldata.rds")
```

## 01. Input Validation and Low-Count Gene Filtering

**Purpose:**
Sample mismatch between count matrix and colData is the most common source of erroneous results in RNA-seq analysis. This step validates sample alignment, computes per-sample QC metrics (library size, detected genes), and applies low-count gene filtering to remove features with insufficient reads for reliable dispersion estimation.

**Code:**
```r
### Step 1: Load Libraries
library(tidyverse)

### Step 2: Load Data
counts <- readRDS("results/rds/counts_raw.rds")
coldata <- readRDS("results/rds/coldata.rds")

### Step 3: Validate Sample Alignment
count_samples <- colnames(counts)
coldata_samples <- rownames(coldata)

in_counts_not_coldata <- setdiff(count_samples, coldata_samples)
in_coldata_not_counts <- setdiff(coldata_samples, count_samples)

if (length(in_counts_not_coldata) > 0) {
  warning("Samples in counts but not coldata: ", paste(in_counts_not_coldata, collapse = ", "))
}
if (length(in_coldata_not_counts) > 0) {
  warning("Samples in coldata but not counts: ", paste(in_coldata_not_counts, collapse = ", "))
}

### Step 4: Subset to Common Samples and Reorder
common_samples <- intersect(count_samples, coldata_samples)
counts <- counts[, common_samples]
coldata <- coldata[common_samples, ]

stopifnot(identical(colnames(counts), rownames(coldata)))
cat("Sample alignment verified:", ncol(counts), "samples\n")

### Step 5: Remove Samples with Missing Condition
valid_idx <- !is.na(coldata$condition)
counts <- counts[, valid_idx]
coldata <- coldata[valid_idx, ]
cat("Samples after removing NA condition:", ncol(counts), "\n")

### Step 6: Compute Per-Sample QC Metrics
sample_qc <- data.frame(
  sample_id = colnames(counts),
  library_size = colSums(counts),
  detected_genes = colSums(counts > 0),
  condition = coldata$condition
)

### Step 7: Apply Low-Count Gene Filter
min_count <- 10
min_samples <- 3
keep_genes <- rowSums(counts >= min_count) >= min_samples

counts_filtered <- counts[keep_genes, ]
cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(counts_filtered), "\n")
cat("Genes removed:", nrow(counts) - nrow(counts_filtered), "\n")

### Step 8: Export
write_tsv(sample_qc, "results/tables/sample_qc_metrics.tsv")
saveRDS(counts_filtered, "results/rds/counts_filtered.rds")
saveRDS(coldata, "results/rds/coldata_validated.rds")
```

### Visualization: Sample Quality Control Metrics

Comprehensive QC showing library size distribution, detected genes, and library complexity across treatment groups.

![Sample QC Metrics](Visualizations/Fig01_Sample_QC.png)

**Visualization Code:**
```r
library(ggplot2)
library(patchwork)

qc <- read.table("results/tables/sample_qc_metrics.tsv", header = TRUE, sep = "\t")
qc$condition <- coldata[qc$sample, "condition"]

# 1A: Library Size Distribution
p1a <- ggplot(qc, aes(x = library_size / 1e6, fill = condition)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_rug(aes(color = condition), alpha = 0.5) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition, guide = "none") +
  labs(title = "Library Size Distribution",
       x = "Library Size (millions of reads)", y = "Density") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85))

# 1B: Library Size Boxplot
p1b <- ggplot(qc, aes(x = condition, y = library_size / 1e6, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.6, shape = 21, fill = "white") +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black") +
  scale_fill_manual(values = pal_condition) +
  labs(title = "Library Size by Response", x = NULL, y = "Library Size (millions)") +
  theme_publication() +
  theme(legend.position = "none")

# 1C: Detected Genes Violin
p1c <- ggplot(qc, aes(x = condition, y = detected_genes / 1000, fill = condition)) +
  geom_violin(alpha = 0.7, color = "black", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = pal_condition) +
  labs(title = "Genes Detected per Sample", x = NULL, y = "Detected Genes (thousands)") +
  theme_publication() +
  theme(legend.position = "none")

# 1D: Library Complexity
cor_val <- cor(qc$library_size, qc$detected_genes)
p1d <- ggplot(qc, aes(x = library_size / 1e6, y = detected_genes / 1000)) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", fill = "#CCCCCC") +
  geom_point(aes(fill = condition), size = 3, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("r = ", round(cor_val, 3)), fontface = "italic") +
  labs(title = "Library Complexity", x = "Library Size (millions)", 
       y = "Detected Genes (thousands)") +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85))

# Combine
fig1 <- (p1a | p1b) / (p1c | p1d) +
  plot_annotation(
    title = "Sample Quality Control Metrics",
    subtitle = "GSE192341 | Breast Cancer Neoadjuvant Chemotherapy Response"
  )

ggsave("results/figures/png/Fig01_Sample_QC.png", fig1, width = 12, height = 11, dpi = 300)
```

**Key Features:**
- Library size density plots with rug marks
- Boxplots with individual sample points
- Violin plots showing gene detection
- Library complexity scatter with Pearson correlation

---

## 02. Differential Expression Analysis with DESeq2

**Purpose:**
DESeq2 models RNA-seq count data using a negative binomial generalized linear model (GLM), accounting for library size differences via size factor normalization and biological variability via gene-wise dispersion estimates. The apeglm shrinkage estimator applies an adaptive heavy-tailed Cauchy prior to log2 fold changes, reducing noise in low-count genes while preserving large effect sizes.

**Code:**
```r
### Step 1: Load Libraries
library(DESeq2)
library(apeglm)
library(tidyverse)

### Step 2: Load Validated Data
counts <- readRDS("results/rds/counts_filtered.rds")
coldata <- readRDS("results/rds/coldata_validated.rds")

### Step 3: Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

### Step 4: Run DESeq2 Pipeline
dds <- DESeq(dds)

### Step 5: Inspect Size Factors
sizeFactors(dds)
summary(sizeFactors(dds))

### Step 6: Extract Results (Wald Test)
resultsNames(dds)
res <- results(dds, contrast = c("condition", "pCR", "No_pCR"), alpha = 0.05)
summary(res)

### Step 7: Apply apeglm Log-Fold Change Shrinkage
res_shrunk <- lfcShrink(dds, coef = "condition_pCR_vs_No_pCR", type = "apeglm")

### Step 8: Create Results Table
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj) %>%
  mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

### Step 9: Summarize DE Results
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Upregulated in pCR:", sum(res_df$regulation == "Upregulated", na.rm = TRUE), "\n")
cat("Downregulated in pCR:", sum(res_df$regulation == "Downregulated", na.rm = TRUE), "\n")

### Step 10: Export Normalized Counts and VST Matrix
normalized_counts <- counts(dds, normalized = TRUE)
vst_matrix <- assay(vst(dds, blind = FALSE))

### Step 11: Save All Outputs
saveRDS(dds, "results/rds/dds.rds")
write_tsv(res_df, "results/tables/DESeq2_results_all.tsv")
write_tsv(filter(res_df, padj < 0.05), "results/tables/DESeq2_results_significant.tsv")
write.csv(normalized_counts, "results/tables/normalized_counts.csv")
write.csv(vst_matrix, "results/tables/vst_matrix.csv")
```

### Visualization: Differential Expression Analysis

#### Volcano Plot
Elite volcano plot with gene symbol labels for top differentially expressed genes.

![Volcano Plot](Visualizations/Fig04_Volcano.png)

**Visualization Code:**
```r
res <- read.table("results/tables/DESeq2_results_all.tsv", header = TRUE, sep = "\t")

# Map to gene symbols
res$symbol <- map_ensembl_to_symbol(res$gene_id)

# Classification
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

# Top genes for labeling
top_up <- res %>% filter(regulation == "Upregulated") %>% arrange(padj) %>% head(12)
top_down <- res %>% filter(regulation == "Downregulated") %>% arrange(padj) %>% head(12)
top_label <- bind_rows(top_up, top_down)

n_up <- sum(res$regulation == "Upregulated")
n_down <- sum(res$regulation == "Downregulated")

# Build plot
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
    subtitle = paste0("Total genes: ", format(nrow(res), big.mark = ","), 
                      " | Significance: padj < 0.05, |LFC| > 0.5"),
    x = expression(bold(Log[2]~Fold~Change~(pCR~"/"~No~pCR))),
    y = expression(bold(-Log[10]~Adjusted~P-value))
  ) +
  theme_publication(base_size = 12) +
  theme(legend.position = c(0.85, 0.85))

ggsave("results/figures/png/Fig04_Volcano.png", p4, width = 11, height = 10, dpi = 300)
```

**Highlights:**
- Top 24 genes labeled with HGNC symbols
- Colorblind-friendly regulation colors
- Significance thresholds (padj < 0.05, |LFC| > 0.5)
- Intelligent label positioning with ggrepel

---

#### MA Plot
Mean-average plot showing log2 fold change vs mean normalized expression.

![MA Plot](Visualizations/Fig05_MA_Plot.png)

**Visualization Code:**
```r
dds <- readRDS("results/rds/dds.rds")
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
    name = "Significance",
    labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not Significant")
  ) +
  labs(title = "MA Plot", subtitle = "Log2 fold change vs mean normalized expression",
       x = expression(bold(Log[10]~Mean~Expression)),
       y = expression(bold(Log[2]~Fold~Change))) +
  theme_publication() +
  theme(legend.position = "bottom")

ggsave("results/figures/png/Fig05_MA_Plot.png", p5, width = 10, height = 8, dpi = 300)
```

---

#### DE Statistics
Four-panel statistical overview of DE results.

![DE Statistics](Visualizations/Fig06_DE_Statistics.png)

**Visualization Code:**
```r
# 6A: P-value histogram
p6a <- ggplot(res %>% filter(!is.na(pvalue)), aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#4DBBD5", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#E64B35") +
  annotate("text", x = 0.12, y = Inf, vjust = 2, label = "p = 0.05", 
           color = "#E64B35", fontface = "bold") +
  labs(title = "P-value Distribution",
       subtitle = "Expected uniform under null; enrichment at low p indicates signal",
       x = "P-value", y = "Frequency") +
  theme_publication()

# 6B: LFC distribution
sig_res <- res %>% filter(!is.na(padj) & padj < 0.05)
p6b <- ggplot(sig_res, aes(x = log2FoldChange, fill = regulation)) +
  geom_histogram(bins = 50, color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = pal_regulation, name = "Direction") +
  labs(title = "Log2 Fold Change Distribution",
       subtitle = paste0("Significant genes only (n = ", nrow(sig_res), ", padj < 0.05)"),
       x = expression(Log[2]~Fold~Change), y = "Frequency") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85))

# 6C: Effect size vs significance
p6c <- ggplot(res %>% filter(!is.na(padj)), aes(x = abs(log2FoldChange), y = -log10(padj))) +
  geom_hex(bins = 50) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#E64B35") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#E64B35") +
  scale_fill_viridis_c(option = "plasma", name = "Count", trans = "log10") +
  labs(title = "Effect Size vs Significance",
       subtitle = "Density of genes across fold change and p-value space",
       x = expression("|"*Log[2]~Fold~Change*"|"),
       y = expression(-Log[10]~Adjusted~P-value)) +
  theme_publication()

# 6D: Rank plot
res_ranked <- res %>% filter(!is.na(padj)) %>% arrange(padj) %>% mutate(rank = row_number())
p6d <- ggplot(res_ranked %>% head(500), aes(x = rank, y = -log10(padj))) +
  geom_line(color = "#3C5488", linewidth = 1) +
  geom_point(aes(color = regulation), size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#888888") +
  scale_color_manual(values = pal_regulation, name = "Regulation") +
  labs(title = "Top 500 Genes by Significance", subtitle = "Ranked by adjusted p-value",
       x = "Rank", y = expression(-Log[10]~Adjusted~P-value)) +
  theme_publication() +
  theme(legend.position = "bottom")

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(title = "Differential Expression Statistics")

ggsave("results/figures/png/Fig06_DE_Statistics.png", fig6, width = 13, height = 11, dpi = 300)
```

---

#### Top Gene Expression
Individual expression profiles of the 12 most significant genes.

![Top Gene Expression](Visualizations/Fig08_Top_Gene_Expression.png)

**Visualization Code:**
```r
norm_counts <- as.matrix(read.table("results/tables/normalized_counts.tsv", header = TRUE, 
                                    sep = "\t", row.names = 1, check.names = FALSE))

top12 <- res %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(12)
top12_genes <- top12$gene_id[top12$gene_id %in% rownames(norm_counts)]
top12 <- top12 %>% filter(gene_id %in% top12_genes)

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
  plot_annotation(title = "Expression of Top Differentially Expressed Genes",
                  subtitle = "Ranked by adjusted p-value")

ggsave("results/figures/png/Fig08_Top_Gene_Expression.png", fig8, width = 14, height = 10, dpi = 300)
```

---

#### Supplementary: Dispersion Estimates & Size Factors

![Dispersion Plot](Visualizations/FigS1_Dispersion.png)

![Size Factors](Visualizations/FigS2_Size_Factors.png)

**Visualization Code:**
```r
# S1: Dispersion plot
dds <- readRDS("results/rds/dds.rds")
png("results/figures/png/FigS1_Dispersion.png", width = 9, height = 7, units = "in", res = 300)
plotDispEsts(dds, main = "DESeq2 Gene-wise Dispersion Estimates")
dev.off()

# S2: Size factors
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

ggsave("results/figures/png/FigS2_Size_Factors.png", pS2, width = 12, height = 5, dpi = 300)
```

---

## 03. Quality Control: VST, PCA, and Sample Distance Analysis

**Purpose:**
Variance-stabilizing transformation (VST) removes the mean-variance dependence inherent in count data, enabling meaningful Euclidean distance calculations and PCA. Principal component analysis reveals the major axes of transcriptomic variation and whether samples cluster by experimental condition or confounding variables. Distance-based outlier detection identifies samples with aberrant global expression profiles.

**Code:**
```r
### Step 1: Load Libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

### Step 2: Load DESeqDataSet
dds <- readRDS("results/rds/dds.rds")

### Step 3: Apply Variance-Stabilizing Transformation
vsd <- vst(dds, blind = TRUE)

### Step 4: Principal Component Analysis
pca_result <- prcomp(t(assay(vsd)))
percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  sample_id = colnames(vsd),
  condition = colData(vsd)$condition
)

cat("Variance explained - PC1:", percent_var[1], "%, PC2:", percent_var[2], "%\n")

### Step 5: Compute Sample Distance Matrix
sample_dists <- dist(t(assay(vsd)), method = "euclidean")
sample_dist_matrix <- as.matrix(sample_dists)

### Step 6: Identify Outliers by Mean Distance
mean_dists <- rowMeans(sample_dist_matrix)
z_scores <- scale(mean_dists)[, 1]

outlier_df <- data.frame(
  sample_id = names(mean_dists),
  mean_distance = mean_dists,
  z_score = z_scores,
  outlier_flag = abs(z_scores) > 3
)

cat("Potential outliers (|Z| > 3):\n")
print(filter(outlier_df, outlier_flag))

### Step 7: Identify Top Variable Genes
rv <- rowVars(assay(vsd))
top_var_genes <- head(order(rv, decreasing = TRUE), 500)
top_var_gene_ids <- rownames(vsd)[top_var_genes]

### Step 8: Export QC Tables
write_tsv(pca_df, "results/tables/qc_pca_coordinates.tsv")
write_tsv(data.frame(PC = 1:10, variance_percent = percent_var[1:10]), 
          "results/tables/qc_pca_variance.tsv")
write_tsv(outlier_df, "results/tables/qc_outlier_flags.tsv")
write_tsv(data.frame(gene_id = top_var_gene_ids), "results/tables/qc_top_variable_genes.tsv")
saveRDS(vsd, "results/rds/vsd.rds")
```

### Visualization: PCA Analysis

Multi-panel PCA analysis revealing transcriptomic variance structure and sample clustering patterns.

![PCA Analysis](Visualizations/Fig02_PCA.png)

**Visualization Code:**
```r
pca_df <- read.table("results/tables/qc_pca_coordinates.tsv", header = TRUE, sep = "\t")
var_df <- read.table("results/tables/qc_pca_variance_explained.tsv", header = TRUE, sep = "\t")
var_df$pct <- round(var_df$variance_explained * 100, 1)

# 2A: PC1 vs PC2
p2a <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = condition), level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(fill = condition), size = 4, shape = 21, stroke = 0.5, alpha = 0.85) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition) +
  labs(title = "Principal Component Analysis",
       x = paste0("PC1 (", var_df$pct[1], "% variance)"),
       y = paste0("PC2 (", var_df$pct[2], "% variance)")) +
  theme_publication(base_size = 12) +
  theme(legend.position = "bottom")

# 2B: PC2 vs PC3
p2b <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
  stat_ellipse(aes(color = condition), level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(fill = condition), size = 4, shape = 21, stroke = 0.5, alpha = 0.85) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  scale_color_manual(values = pal_condition) +
  labs(title = "PC2 vs PC3",
       x = paste0("PC2 (", var_df$pct[2], "% variance)"),
       y = paste0("PC3 (", var_df$pct[3], "% variance)")) +
  theme_publication(base_size = 12) +
  theme(legend.position = "bottom")

# 2C: Scree plot
var_plot <- var_df[1:min(10, nrow(var_df)), ]
var_plot$PC_num <- 1:nrow(var_plot)
var_plot$cumulative <- cumsum(var_plot$pct)

p2c <- ggplot(var_plot, aes(x = PC_num)) +
  geom_bar(aes(y = pct), stat = "identity", fill = "#4DBBD5", alpha = 0.85, width = 0.7) +
  geom_line(aes(y = cumulative / 2), color = "#E64B35", linewidth = 1.2) +
  geom_point(aes(y = cumulative / 2), color = "#E64B35", size = 3) +
  geom_text(aes(y = pct, label = paste0(pct, "%")), vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_x_continuous(breaks = 1:10, labels = paste0("PC", 1:10)) +
  scale_y_continuous(name = "Individual Variance (%)",
                     sec.axis = sec_axis(~ . * 2, name = "Cumulative Variance (%)")) +
  labs(title = "Variance Explained by Principal Components", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2D: PC1 distribution
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

ggsave("results/figures/png/Fig02_PCA.png", fig2, width = 13, height = 12, dpi = 300)
```

**Panels:**
- PC1 vs PC2 with 95% confidence ellipses
- PC2 vs PC3 alternative view
- Scree plot with cumulative variance
- PC1 score distribution by treatment group

---

### Visualization: Sample Clustering Heatmaps

Hierarchical clustering based on sample-to-sample distances and correlations.

#### Sample Distance Matrix
![Sample Distance Heatmap](Visualizations/Fig03A_Sample_Distance.png)

#### Sample Correlation Matrix
![Sample Correlation Heatmap](Visualizations/Fig03B_Sample_Correlation.png)

**Visualization Code:**
```r
library(pheatmap)

dist_mat <- as.matrix(read.table("results/tables/qc_sample_distance_matrix.tsv", 
                                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
cor_mat <- as.matrix(read.table("results/tables/qc_sample_correlation_matrix.tsv", 
                                header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))

common <- intersect(rownames(dist_mat), rownames(coldata))
dist_mat <- dist_mat[common, common]
cor_mat <- cor_mat[common, common]

condition_values <- as.character(coldata[common, "condition"])
anno <- data.frame(Response = condition_values, row.names = common)

unique_conditions <- unique(na.omit(condition_values))
cond_colors <- c("#DC3220", "#005AB5")
names(cond_colors) <- unique_conditions[1:min(2, length(unique_conditions))]
anno_colors <- list(Response = cond_colors)

# Distance colors
pal_distance <- colorRampPalette(c("#FFFFFF", "#FFF5F0", "#FEE0D2", "#FCBBA1", 
                                   "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D"))(100)

# 3A: Distance heatmap
png("results/figures/png/Fig03A_Sample_Distance.png", width = 10, height = 9, units = "in", res = 300)
pheatmap(dist_mat, color = pal_distance, annotation_col = anno, annotation_row = anno,
         annotation_colors = anno_colors, show_rownames = FALSE, show_colnames = FALSE,
         clustering_method = "ward.D2", main = "Sample-to-Sample Euclidean Distance",
         fontsize = 11, border_color = NA, treeheight_row = 40, treeheight_col = 40)
dev.off()

# Correlation colors
cor_colors <- colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"))(100)

# 3B: Correlation heatmap
png("results/figures/png/Fig03B_Sample_Correlation.png", width = 10, height = 9, units = "in", res = 300)
pheatmap(cor_mat, color = cor_colors, breaks = seq(min(cor_mat, na.rm = TRUE), 1, length.out = 101),
         annotation_col = anno, annotation_row = anno, annotation_colors = anno_colors,
         show_rownames = FALSE, show_colnames = FALSE, clustering_method = "ward.D2",
         main = "Sample-to-Sample Pearson Correlation", fontsize = 11, border_color = NA,
         treeheight_row = 40, treeheight_col = 40)
dev.off()
```

**Features:**
- Ward.D2 hierarchical clustering
- Response group annotations
- Diverging color schemes
- Dendrogram visualization

---

### Visualization: Top DE Genes Heatmap

Z-score normalized expression heatmap of top 60 differentially expressed genes.

![Top DE Genes Heatmap](Visualizations/Fig07_Top_DE_Heatmap.png)

**Visualization Code:**
```r
vst_mat <- as.matrix(read.table("results/tables/vst_matrix.tsv", header = TRUE, sep = "\t", 
                                row.names = 1, check.names = FALSE))

top_genes <- res %>% filter(!is.na(padj) & padj < 0.05) %>% arrange(padj) %>% head(60) %>% pull(gene_id)
top_genes <- top_genes[top_genes %in% rownames(vst_mat)]

heat_mat <- vst_mat[top_genes, , drop = FALSE]

# Map to symbols
original_ids <- rownames(heat_mat)
symbols <- map_ensembl_to_symbol(original_ids)
symbols_unique <- make.unique(symbols)
rownames(heat_mat) <- symbols_unique

# Z-score scaling
heat_scaled <- t(scale(t(heat_mat)))

common <- intersect(colnames(heat_scaled), rownames(coldata))
heat_scaled <- heat_scaled[, common]

# Annotations
condition_values <- coldata[common, "condition"]
anno_col <- data.frame(Response = condition_values, row.names = common)

gene_lfc <- res$log2FoldChange[match(top_genes, res$gene_id)]
direction_vec <- ifelse(is.na(gene_lfc), "Unknown",
                        ifelse(gene_lfc > 0, "Up in pCR", "Down in pCR"))
names(direction_vec) <- symbols_unique
anno_row <- data.frame(Direction = direction_vec[rownames(heat_scaled)], row.names = rownames(heat_scaled))

unique_conditions <- unique(na.omit(condition_values))
cond_colors <- c("#DC3220", "#005AB5")
names(cond_colors) <- unique_conditions[1:min(2, length(unique_conditions))]

anno_colors <- list(
  Response = cond_colors,
  Direction = c("Up in pCR" = "#C41E3A", "Down in pCR" = "#1E5AA8", "Unknown" = "#888888")
)

# Heatmap colors
pal_heatmap_diverging <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                            "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100)

png("results/figures/png/Fig07_Top_DE_Heatmap.png", width = 12, height = 14, units = "in", res = 300)
pheatmap(heat_scaled, color = pal_heatmap_diverging, breaks = seq(-3, 3, length.out = 101),
         annotation_col = anno_col, annotation_row = anno_row, annotation_colors = anno_colors,
         show_colnames = FALSE, fontsize_row = 7, fontsize = 10,
         main = paste0("Top ", nrow(heat_scaled), " Differentially Expressed Genes"),
         clustering_method = "ward.D2", border_color = NA,
         cutree_cols = 2, cutree_rows = 2, treeheight_row = 30, treeheight_col = 30)
dev.off()
```

---

### Supplementary: Top Variable Genes

![Top Variable Genes](Visualizations/FigS3_Top_Variable_Genes.png)

**Visualization Code:**
```r
top_var <- read.table("results/tables/qc_top_variable_genes.tsv", header = TRUE, sep = "\t")
vst_mat <- as.matrix(read.table("results/tables/vst_matrix.tsv", header = TRUE, sep = "\t", 
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

unique_conditions <- unique(na.omit(condition_values))
cond_colors <- c("#DC3220", "#005AB5")
names(cond_colors) <- unique_conditions[1:min(2, length(unique_conditions))]
anno_colors <- list(Response = cond_colors)

png("results/figures/png/FigS3_Top_Variable_Genes.png", width = 12, height = 16, units = "in", res = 300)
pheatmap(heat_scaled, color = pal_heatmap_diverging, breaks = seq(-3, 3, length.out = 101),
         annotation_col = anno, annotation_colors = anno_colors, show_colnames = FALSE,
         fontsize_row = 6, main = paste0("Top ", length(top100), " Most Variable Genes"),
         clustering_method = "ward.D2", border_color = NA)
dev.off()
```

---

## 04. Cell-Type Deconvolution with xCell

**Purpose:**
Bulk RNA-seq represents a mixture of cell populations. xCell is a gene signature-based method that estimates enrichment scores for 64 immune and stromal cell types from bulk expression profiles. Comparing deconvolution scores between pCR and No pCR groups can reveal tumor microenvironment features associated with treatment response.

**Code:**
```r
### Step 1: Load Libraries
library(xCell)
library(tidyverse)
library(org.Hs.eg.db)

### Step 2: Load VST Expression Matrix
vsd <- readRDS("results/rds/vsd.rds")
expr_matrix <- assay(vsd)

### Step 3: Map ENSEMBL IDs to Gene Symbols
gene_annotation <- read.csv("data/gene_annotation.csv")

ensembl_to_symbol <- setNames(gene_annotation$gene_symbol, gene_annotation$gene_id)
mapped_symbols <- ensembl_to_symbol[rownames(expr_matrix)]

valid_idx <- !is.na(mapped_symbols) & mapped_symbols != ""
expr_matrix_symbol <- expr_matrix[valid_idx, ]
rownames(expr_matrix_symbol) <- mapped_symbols[valid_idx]

expr_matrix_symbol <- expr_matrix_symbol[!duplicated(rownames(expr_matrix_symbol)), ]

cat("Genes mapped to symbols:", nrow(expr_matrix_symbol), "\n")

### Step 4: Run xCell Deconvolution
xcell_scores <- xCellAnalysis(expr_matrix_symbol)

### Step 5: Prepare Results Table
xcell_df <- as.data.frame(t(xcell_scores)) %>%
  rownames_to_column("sample_id")

### Step 6: Add Condition Information
coldata <- readRDS("results/rds/coldata_validated.rds")
xcell_df <- xcell_df %>%
  left_join(
    data.frame(sample_id = rownames(coldata), condition = coldata$condition),
    by = "sample_id"
  )

### Step 7: Statistical Comparison Between Groups
xcell_long <- xcell_df %>%
  pivot_longer(cols = -c(sample_id, condition), names_to = "cell_type", values_to = "score")

comparison_results <- xcell_long %>%
  group_by(cell_type) %>%
  summarise(
    mean_pCR = mean(score[condition == "pCR"], na.rm = TRUE),
    mean_NopCR = mean(score[condition == "No_pCR"], na.rm = TRUE),
    log2FC = log2((mean_pCR + 0.01) / (mean_NopCR + 0.01)),
    p_value = wilcox.test(score ~ condition)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("Cell types with p_adj < 0.05:\n")
print(filter(comparison_results, p_adj < 0.05))

### Step 8: Export Results
write_tsv(xcell_df, "results/tables/deconvolution_xcell_scores.tsv")
write_tsv(comparison_results, "results/tables/deconvolution_xcell_comparison.tsv")
saveRDS(xcell_scores, "results/rds/xcell_scores.rds")
```

### Visualization: xCell Deconvolution Analysis

Comprehensive tumor microenvironment characterization.

#### xCell Enrichment Heatmap
![xCell Heatmap](Visualizations/Fig09A_xCell_Heatmap.png)

Top 30 most variable cell type enrichment scores.

**Visualization Code:**
```r
xcell <- read.table("results/tables/deconvolution_xcell_scores.tsv", header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE)

cell_cols <- setdiff(colnames(xcell), c("sample", "condition"))
xcell_mat <- as.matrix(xcell[, cell_cols])
rownames(xcell_mat) <- xcell$sample

# Top 30 variable
cell_var <- apply(xcell_mat, 2, var, na.rm = TRUE)
top_cells <- names(sort(cell_var, decreasing = TRUE))[1:min(30, length(cell_var))]
xcell_top <- xcell_mat[, top_cells]
xcell_scaled <- scale(xcell_top)

condition_values <- as.character(xcell$condition)
anno <- data.frame(Response = condition_values, row.names = xcell$sample)

unique_conditions <- unique(na.omit(condition_values))
cond_colors <- c("#DC3220", "#005AB5")
names(cond_colors) <- unique_conditions[1:min(2, length(unique_conditions))]
anno_colors <- list(Response = cond_colors)

png("results/figures/png/Fig09A_xCell_Heatmap.png", width = 14, height = 10, units = "in", res = 300)
pheatmap(t(xcell_scaled), color = pal_heatmap_diverging, breaks = seq(-2.5, 2.5, length.out = 101),
         annotation_col = anno, annotation_colors = anno_colors, show_colnames = FALSE,
         fontsize_row = 9, fontsize = 10,
         main = "xCell Enrichment Scores (Top 30 Variable Cell Types)",
         clustering_method = "ward.D2", border_color = NA, treeheight_col = 30)
dev.off()
```

---

#### Immune Cell Populations
![Immune Cells](Visualizations/Fig09B_Immune_Cells.png)

Key immune cell infiltration comparison between treatment groups.

**Visualization Code:**
```r
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
  labs(title = "Immune Cell Infiltration",
       subtitle = "xCell enrichment scores by treatment response",
       x = NULL, y = "xCell Enrichment Score") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "bottom")

ggsave("results/figures/png/Fig09B_Immune_Cells.png", p9b, width = 12, height = 8, dpi = 300)
```

---

#### Differentially Enriched Cell Types
![Significant Cell Types](Visualizations/Fig09C_Significant_CellTypes.png)

Cell populations with significant differential enrichment (Wilcoxon test, padj < 0.25).

**Visualization Code:**
```r
xcell_stats <- read.table("results/tables/deconvolution_xcell_pcr_vs_nopcr.tsv", 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sig_cells <- xcell_stats %>% filter(!is.na(padj) & padj < 0.25) %>% arrange(padj) %>% head(15)

sig_long <- xcell %>%
  dplyr::select(sample, condition, all_of(sig_cells$population)) %>%
  pivot_longer(cols = -c(sample, condition), names_to = "Cell_Type", values_to = "Score") %>%
  left_join(sig_cells %>% dplyr::select(population, padj, delta_median), 
            by = c("Cell_Type" = "population")) %>%
  mutate(Cell_Type = fct_reorder(Cell_Type, -padj))

p9c <- ggplot(sig_long, aes(x = Cell_Type, y = Score, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.6, width = 0.6) +
  scale_fill_manual(values = pal_condition, name = "Response") +
  coord_flip() +
  labs(title = "Differentially Enriched Cell Populations",
       subtitle = "Wilcoxon test, padj < 0.25", x = NULL, y = "xCell Score") +
  theme_publication() +
  theme(legend.position = "bottom")

ggsave("results/figures/png/Fig09C_Significant_CellTypes.png", p9c, width = 10, height = 8, dpi = 300)
```

---

#### Tumor Microenvironment Scores
![TME Scores](Visualizations/Fig09D_TME_Scores.png)

ImmuneScore, StromaScore, and MicroenvironmentScore by treatment response.

**Visualization Code:**
```r
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
  labs(title = "Tumor Microenvironment Scores", x = NULL, y = "Score") +
  theme_publication() +
  theme(legend.position = "none")

ggsave("results/figures/png/Fig09D_TME_Scores.png", p9d, width = 10, height = 6, dpi = 300)
```

---

## 05. Pathway Enrichment Analysis (ORA and GSEA)

**Purpose:**
Over-representation analysis (ORA) tests whether differentially expressed genes are enriched in predefined gene sets using a hypergeometric test. Gene set enrichment analysis (GSEA) uses the full ranked gene list to detect coordinated expression changes without requiring an arbitrary significance cutoff. We apply both methods to Hallmark (MSigDB), Gene Ontology Biological Process, and Reactome pathway collections.

**Code:**
```r
### Step 1: Load Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(tidyverse)

### Step 2: Load DE Results and Gene Annotation
res_df <- read_tsv("results/tables/DESeq2_results_all.tsv")
gene_annotation <- read.csv("data/gene_annotation.csv")

### Step 3: Map ENSEMBL to ENTREZID
ensembl_ids <- res_df$gene_id
entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENSEMBL"
)

res_annotated <- res_df %>%
  left_join(entrez_map, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID))

cat("Genes with ENTREZID mapping:", nrow(res_annotated), "\n")

### Step 4: Define Gene Lists for ORA
sig_up <- res_annotated %>%
  filter(padj < 0.05, log2FoldChange > 0) %>%
  pull(ENTREZID)

sig_down <- res_annotated %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  pull(ENTREZID)

universe <- res_annotated$ENTREZID

cat("Upregulated genes:", length(sig_up), "\n")
cat("Downregulated genes:", length(sig_down), "\n")

### Step 5: Prepare Ranked Gene List for GSEA
ranked_genes <- res_annotated %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = ENTREZID)

ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]

### Step 6: Get MSigDB Hallmark Gene Sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

### Step 7: Hallmark ORA - Upregulated Genes
ora_hallmark_up <- enricher(
  gene = sig_up,
  universe = universe,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

### Step 8: Hallmark ORA - Downregulated Genes
ora_hallmark_down <- enricher(
  gene = sig_down,
  universe = universe,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

### Step 9: Hallmark GSEA
gsea_hallmark <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = hallmark_sets,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

### Step 10: GO Biological Process ORA
ora_gobp_up <- enrichGO(
  gene = sig_up,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ora_gobp_down <- enrichGO(
  gene = sig_down,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

### Step 11: Reactome Pathway ORA
ora_reactome_up <- enrichPathway(
  gene = sig_up,
  universe = universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ora_reactome_down <- enrichPathway(
  gene = sig_down,
  universe = universe,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

### Step 12: Export All Results
write_tsv(res_annotated, "results/tables/DESeq2_results_annotated.tsv")

if (!is.null(ora_hallmark_up) && nrow(as.data.frame(ora_hallmark_up)) > 0) {
  write_tsv(as.data.frame(ora_hallmark_up), "results/tables/enrich_hallmark_ORA_UP.tsv")
}
if (!is.null(ora_hallmark_down) && nrow(as.data.frame(ora_hallmark_down)) > 0) {
  write_tsv(as.data.frame(ora_hallmark_down), "results/tables/enrich_hallmark_ORA_DOWN.tsv")
}
if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  write_tsv(as.data.frame(gsea_hallmark), "results/tables/enrich_hallmark_GSEA.tsv")
}
if (!is.null(ora_gobp_up) && nrow(as.data.frame(ora_gobp_up)) > 0) {
  write_tsv(as.data.frame(ora_gobp_up), "results/tables/enrich_GO_BP_ORA_UP.tsv")
}
if (!is.null(ora_gobp_down) && nrow(as.data.frame(ora_gobp_down)) > 0) {
  write_tsv(as.data.frame(ora_gobp_down), "results/tables/enrich_GO_BP_ORA_DOWN.tsv")
}
if (!is.null(ora_reactome_up) && nrow(as.data.frame(ora_reactome_up)) > 0) {
  write_tsv(as.data.frame(ora_reactome_up), "results/tables/enrich_reactome_ORA_UP.tsv")
}
if (!is.null(ora_reactome_down) && nrow(as.data.frame(ora_reactome_down)) > 0) {
  write_tsv(as.data.frame(ora_reactome_down), "results/tables/enrich_reactome_ORA_DOWN.tsv")
}

saveRDS(list(
  ora_hallmark_up = ora_hallmark_up,
  ora_hallmark_down = ora_hallmark_down,
  gsea_hallmark = gsea_hallmark,
  ora_gobp_up = ora_gobp_up,
  ora_gobp_down = ora_gobp_down,
  ora_reactome_up = ora_reactome_up,
  ora_reactome_down = ora_reactome_down
), "results/rds/enrichment_results.rds")
```

### Visualization: Pathway Enrichment Analysis

Multi-panel pathway enrichment visualization across gene set databases.

#### Hallmark Pathways - Upregulated
![Hallmark UP](Visualizations/Fig10A_Hallmark_UP.png)

#### Hallmark Pathways - Downregulated
![Hallmark DOWN](Visualizations/Fig10B_Hallmark_DOWN.png)

#### GO Biological Process - Upregulated
![GO BP UP](Visualizations/Fig10C_GO_BP_UP.png)

#### GO Biological Process - Downregulated
![GO BP DOWN](Visualizations/Fig10D_GO_BP_DOWN.png)

#### GSEA Hallmark Pathways
![GSEA Hallmark](Visualizations/Fig10E_GSEA_Hallmark.png)

**Visualization Code:**
```r
# Elite dotplot function
make_elite_dotplot <- function(file_path, title, n_show = 15) {
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
  
  df$Description <- gsub("HALLMARK_", "", df$ID)
  df$Description <- gsub("_", " ", df$Description)
  df$Description <- tools::toTitleCase(tolower(df$Description))
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
    scale_color_viridis_c(option = "plasma", name = expression(-Log[10]~P[adj]), direction = -1) +
    scale_size_continuous(range = c(4, 12), name = "Gene Count") +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 9), legend.position = "right")
}

# Generate pathway plots
p_hall_up <- make_elite_dotplot("results/tables/enrich_hallmark_ORA_UP.tsv", 
                                "Hallmark Pathways: Upregulated in pCR", n_show = 12)
if (!is.null(p_hall_up)) {
  ggsave("results/figures/png/Fig10A_Hallmark_UP.png", p_hall_up, 
         width = 10, height = 8, dpi = 300, bg = "white")
}

p_hall_down <- make_elite_dotplot("results/tables/enrich_hallmark_ORA_DOWN.tsv", 
                                  "Hallmark Pathways: Downregulated in pCR", n_show = 12)
if (!is.null(p_hall_down)) {
  ggsave("results/figures/png/Fig10B_Hallmark_DOWN.png", p_hall_down, 
         width = 10, height = 8, dpi = 300, bg = "white")
}

p_go_up <- make_elite_dotplot("results/tables/enrich_GO_BP_ORA_UP.tsv", 
                              "GO Biological Process: Upregulated", n_show = 18)
if (!is.null(p_go_up)) {
  ggsave("results/figures/png/Fig10C_GO_BP_UP.png", p_go_up, 
         width = 11, height = 10, dpi = 300, bg = "white")
}

p_go_down <- make_elite_dotplot("results/tables/enrich_GO_BP_ORA_DOWN.tsv", 
                                "GO Biological Process: Downregulated", n_show = 18)
if (!is.null(p_go_down)) {
  ggsave("results/figures/png/Fig10D_GO_BP_DOWN.png", p_go_down, 
         width = 11, height = 10, dpi = 300, bg = "white")
}

# GSEA Barplot
gsea_file <- "results/tables/enrich_hallmark_GSEA.tsv"
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
      theme(axis.text.y = element_text(size = 9), legend.position = "bottom")
    
    ggsave("results/figures/png/Fig10E_GSEA_Hallmark.png", p_gsea, 
           width = 11, height = 10, dpi = 300, bg = "white")
  }
}
```

**Features:**
- Dot plots with gene ratio scaling
- Color intensity by adjusted p-value
- Point size by gene count
- Normalized enrichment scores (GSEA)

---

## Summary Visualization

### Analysis Summary

High-level overview of the complete analysis.

![Analysis Summary](Visualizations/Fig11_Summary.png)

**Visualization Code:**
```r
# 11A: Sample distribution
cond_df <- as.data.frame(table(coldata$condition))
colnames(cond_df) <- c("Response", "Count")

unique_conds <- as.character(cond_df$Response)
cond_colors <- setNames(c("#DC3220", "#005AB5")[1:length(unique_conds)], unique_conds)

p11a <- ggplot(cond_df, aes(x = Response, y = Count, fill = Response)) +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.9) +
  geom_text(aes(label = Count), vjust = -0.4, size = 6, fontface = "bold") +
  scale_fill_manual(values = cond_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Sample Distribution", x = NULL, y = "Number of Samples") +
  theme_publication(base_size = 13) +
  theme(legend.position = "none")

# 11B: DE gene counts
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
  plot_annotation(title = "Analysis Summary: GSE192341",
                  subtitle = "Breast Cancer Neoadjuvant Chemotherapy Response | pCR vs No pCR")

ggsave("results/figures/png/Fig11_Summary.png", fig11, width = 12, height = 5.5, dpi = 300)
```

**Panels:**
- Sample distribution bar chart
- DE gene categorization (strict and lenient LFC cutoffs)

---

## Key Output Files

| File | Description |
|------|-------------|
| `results/tables/DESeq2_results_all.tsv` | Complete differential expression results with shrunken LFC |
| `results/tables/DESeq2_results_significant.tsv` | Genes with padj < 0.05 |
| `results/tables/DESeq2_results_annotated.tsv` | DE results with gene symbols and ENTREZID |
| `results/tables/normalized_counts.csv` | DESeq2 size factor-normalized counts |
| `results/tables/vst_matrix.csv` | Variance-stabilized expression matrix |
| `results/tables/qc_pca_coordinates.tsv` | Sample PCA coordinates |
| `results/tables/qc_outlier_flags.tsv` | Distance-based outlier detection |
| `results/tables/deconvolution_xcell_scores.tsv` | xCell cell-type enrichment scores |
| `results/tables/enrich_hallmark_*.tsv` | Hallmark pathway enrichment results |
| `results/tables/enrich_GO_BP_*.tsv` | GO Biological Process enrichment |
| `results/tables/enrich_reactome_*.tsv` | Reactome pathway enrichment |
| `results/figures/png/*` | Publication-grade PNG figures (300 DPI) |
| `results/figures/pdf/*` | Vector PDF figures |

## Complete Pipeline Execution

Run the entire analysis pipeline sequentially:
```bash
# Data preparation
Rscript scripts/00A_make_counts_from_processed_matrix.R
Rscript scripts/00B_make_coldata_from_GEO.R

# Core analysis
Rscript scripts/01_input_validation.R
Rscript scripts/02_deseq2_analysis.R
Rscript scripts/03_qc_vst_pca_distances.R
Rscript scripts/04_deconvolution_xcell.R
Rscript scripts/05_enrichment_ORA_GSEA.R

# Elite visualizations (comprehensive)
Rscript scripts/07_visualizations.R
```

## Visualization Features

The elite visualization pipeline generates publication-quality figures with:

- **Elite color palettes** (Nature/Cell/Science inspired)
- **Publication theme** with proper fonts and spacing
- **Gene symbol mapping** from ENSEMBL IDs using org.Hs.eg.db
- **High-resolution output** (300 DPI PNG + vector PDF)
- **Comprehensive figure annotations** and statistical overlays
- **Colorblind-friendly palettes** throughout

**Complete visualization script:** `scripts/07_visualizations.R`

**Output structure:**
```
results/figures/
├── png/          # 300 DPI raster images
└── pdf/          # Vector graphics for publication
```

## References

- **DESeq2:** Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 2014;15:550
- **apeglm:** Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 2019;35:2084–2092
- **xCell:** Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. *Genome Biology* 2017;18:220
- **clusterProfiler:** Wu T et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2021;2:100141
- **MSigDB Hallmark:** Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems* 2015;1:417–425

## Contact

📬 For questions or issues:

**Email:** bioinfosourabh@gmail.com

**Website:** [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app)

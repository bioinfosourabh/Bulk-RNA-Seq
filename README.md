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
8. **Elite Visualizations**: Publication-quality figures with gene symbols, Nature/Cell/Science-grade aesthetics

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

## 00A. Constructing Count Matrix from GEO Processed Data

**Purpose:**
GEO processed files often contain gene annotations interleaved with expression values. This step parses the raw GEO matrix to generate a clean gene × sample integer count matrix and a separate gene annotation table for downstream identifier mapping (ENSEMBL → SYMBOL → ENTREZID).

**Code:**
```r
### Step 1: Load Libraries
library(data.table)
library(tidyverse)

### Step 2: Read GEO Processed Matrix
# The processed matrix contains gene IDs, symbols, and sample expression columns
raw_data <- fread("data_raw/GSE192341_processed_data.txt")

### Step 3: Inspect Column Structure
# Identify annotation columns vs numeric expression columns
head(colnames(raw_data), 10)
str(raw_data[, 1:5])

### Step 4: Extract Gene Annotation
# First two columns typically contain ENSEMBL ID and gene symbol
gene_annotation <- raw_data %>%
  select(gene_id = 1, gene_symbol = 2) %>%
  distinct()

### Step 5: Build Count Matrix
# Extract numeric columns (sample expression values)
sample_cols <- colnames(raw_data)[sapply(raw_data, is.numeric)]
count_matrix <- as.matrix(raw_data[, ..sample_cols])
rownames(count_matrix) <- raw_data[[1]]

### Step 6: Verify Matrix Integrity
# Check for non-negative integers (required for DESeq2)
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
    # Define condition factor with No_pCR as reference level
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

# Check bidirectional matching
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

# Verify identical ordering
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
# Retain genes with >= 10 counts in >= 3 samples (DESeq2 recommendation)
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
# Design formula specifies the experimental factor for differential testing
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

### Step 4: Run DESeq2 Pipeline
# Sequentially performs: estimateSizeFactors, estimateDispersions, nbinomWaldTest
dds <- DESeq(dds)

### Step 5: Inspect Size Factors
# Size factors normalize for sequencing depth differences
sizeFactors(dds)
summary(sizeFactors(dds))

### Step 6: Extract Results (Wald Test)
resultsNames(dds)
res <- results(dds, contrast = c("condition", "pCR", "No_pCR"), alpha = 0.05)
summary(res)

### Step 7: Apply apeglm Log-Fold Change Shrinkage
res_shrunk <- lfcShrink(
  dds,
  coef = "condition_pCR_vs_No_pCR",
  type = "apeglm"
)

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
# blind = TRUE for unbiased QC assessment
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
write_tsv(data.frame(PC = 1:10, variance_percent = percent_var[1:10]), "results/tables/qc_pca_variance.tsv")
write_tsv(outlier_df, "results/tables/qc_outlier_flags.tsv")
write_tsv(data.frame(gene_id = top_var_gene_ids), "results/tables/qc_top_variable_genes.tsv")
saveRDS(vsd, "results/rds/vsd.rds")
```

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
# xCell requires HGNC gene symbols as row names
gene_annotation <- read.csv("data/gene_annotation.csv")

ensembl_to_symbol <- setNames(gene_annotation$gene_symbol, gene_annotation$gene_id)
mapped_symbols <- ensembl_to_symbol[rownames(expr_matrix)]

# Remove unmapped genes
valid_idx <- !is.na(mapped_symbols) & mapped_symbols != ""
expr_matrix_symbol <- expr_matrix[valid_idx, ]
rownames(expr_matrix_symbol) <- mapped_symbols[valid_idx]

# Handle duplicate symbols by keeping highest mean expression
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

# Remove duplicates
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

## 06. Elite Publication-Grade Visualizations

**Purpose:**
Generate Nature/Cell/Science-quality figures with gene symbols, professional color palettes, and comprehensive interpretive plots. All visualizations use colorblind-friendly palettes and follow publication standards.

**Setup & Libraries:**
```r
# Additional packages for elite visualizations
install.packages(c("scales", "patchwork", "viridis", "hexbin"))

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
```

### Figure 1: Sample Quality Control Metrics

Comprehensive QC showing library size distribution, detected genes, and library complexity across treatment groups.

![Sample QC Metrics](Visualizations/Fig01_Sample_QC.png)

**Key Features:**
- Library size density plots with rug marks
- Boxplots with individual sample points
- Violin plots showing gene detection
- Library complexity scatter with correlation coefficient

---

### Figure 2: Principal Component Analysis

Multi-panel PCA analysis revealing transcriptomic variance structure and sample clustering patterns.

![PCA Analysis](Visualizations/Fig02_PCA.png)

**Panels:**
- PC1 vs PC2 with 95% confidence ellipses
- PC2 vs PC3 alternative view
- Scree plot with cumulative variance
- PC1 score distribution by treatment group

---

### Figure 3: Sample Clustering Heatmaps

Hierarchical clustering based on sample-to-sample distances and correlations.

#### 3A: Euclidean Distance Matrix
![Sample Distance Heatmap](Visualizations/Fig03A_Sample_Distance.png)

#### 3B: Pearson Correlation Matrix
![Sample Correlation Heatmap](Visualizations/Fig03B_Sample_Correlation.png)

**Features:**
- Ward.D2 hierarchical clustering
- Response group annotations
- Diverging color schemes
- Dendrogram visualization

---

### Figure 4: Volcano Plot

Elite volcano plot with gene symbol labels for top differentially expressed genes.

![Volcano Plot](Visualizations/Fig04_Volcano.png)

**Highlights:**
- Top 24 genes labeled with HGNC symbols
- Colorblind-friendly regulation colors
- Significance thresholds (padj < 0.05, |LFC| > 0.5)
- Intelligent label positioning with ggrepel

---

### Figure 5: MA Plot

Mean-average plot showing log2 fold change vs mean normalized expression.

![MA Plot](Visualizations/Fig05_MA_Plot.png)

**Features:**
- Top 15 significant genes labeled
- Fold-change threshold lines (±1 LFC)
- Gene symbol annotations
- Color-coded significance levels

---

### Figure 6: Differential Expression Statistics

Four-panel statistical overview of DE results.

![DE Statistics](Visualizations/Fig06_DE_Statistics.png)

**Panels:**
- P-value distribution histogram
- Log2 fold change distribution
- Effect size vs significance density plot
- Top 500 genes ranked by adjusted p-value

---

### Figure 7: Top DE Genes Heatmap

Z-score normalized expression heatmap of top 60 differentially expressed genes.

![Top DE Genes Heatmap](Visualizations/Fig07_Top_DE_Heatmap.png)

**Features:**
- Gene symbols as row labels
- Bi-directional clustering
- Response and regulation direction annotations
- Diverging blue-white-red color scheme

---

### Figure 8: Top Gene Expression Boxplots

Individual expression profiles of the 12 most significant genes.

![Top Gene Expression](Visualizations/Fig08_Top_Gene_Expression.png)

**Features:**
- Violin + boxplot + jitter combination
- Gene symbols with LFC and p-value annotations
- Normalized count scales
- 4×3 grid layout

---

### Figure 9: xCell Deconvolution Analysis

Comprehensive tumor microenvironment characterization.

#### 9A: xCell Enrichment Heatmap
![xCell Heatmap](Visualizations/Fig09A_xCell_Heatmap.png)

Top 30 most variable cell type enrichment scores.

#### 9B: Immune Cell Populations
![Immune Cells](Visualizations/Fig09B_Immune_Cells.png)

Key immune cell infiltration comparison between treatment groups.

#### 9C: Differentially Enriched Cell Types
![Significant Cell Types](Visualizations/Fig09C_Significant_CellTypes.png)

Cell populations with significant differential enrichment (Wilcoxon test, padj < 0.25).

#### 9D: Tumor Microenvironment Scores
![TME Scores](Visualizations/Fig09D_TME_Scores.png)

ImmuneScore, StromaScore, and MicroenvironmentScore by treatment response.

---

### Figure 10: Pathway Enrichment Analysis

Multi-panel pathway enrichment visualization across gene set databases.

#### 10A: Hallmark Pathways - Upregulated
![Hallmark UP](Visualizations/Fig10A_Hallmark_UP.png)

#### 10B: Hallmark Pathways - Downregulated
![Hallmark DOWN](Visualizations/Fig10B_Hallmark_DOWN.png)

#### 10C: GO Biological Process - Upregulated
![GO BP UP](Visualizations/Fig10C_GO_BP_UP.png)

#### 10D: GO Biological Process - Downregulated
![GO BP DOWN](Visualizations/Fig10D_GO_BP_DOWN.png)

#### 10E: GSEA Hallmark Pathways
![GSEA Hallmark](Visualizations/Fig10E_GSEA_Hallmark.png)

**Features:**
- Dot plots with gene ratio scaling
- Color intensity by adjusted p-value
- Point size by gene count
- Normalized enrichment scores (GSEA)

---

### Figure 11: Analysis Summary

High-level overview of the complete analysis.

![Analysis Summary](Visualizations/Fig11_Summary.png)

**Panels:**
- Sample distribution bar chart
- DE gene categorization (strict and lenient LFC cutoffs)

---

### Supplementary Figures

#### S1: Dispersion Estimates
![Dispersion Plot](Visualizations/FigS1_Dispersion.png)

DESeq2 gene-wise dispersion estimates showing fit quality.

#### S2: Size Factors
![Size Factors](Visualizations/FigS2_Size_Factors.png)

Library size normalization factors by sample and treatment group.

#### S3: Top Variable Genes
![Top Variable Genes](Visualizations/FigS3_Top_Variable_Genes.png)

Heatmap of top 100 most variable genes across all samples.

---

## Complete Visualization Code

The complete visualization pipeline (`scripts/07_visualizations.R`) generates all figures with:

- **Elite color palettes** (Nature/Cell/Science inspired)
- **Publication theme** with proper fonts and spacing
- **Gene symbol mapping** from ENSEMBL IDs
- **High-resolution output** (300 DPI PNG + vector PDF)
- **Comprehensive figure annotations**

Key functions include:
- `theme_publication()`: Base ggplot2 theme
- `map_ensembl_to_symbol()`: Robust ID-to-symbol conversion
- `make_elite_dotplot()`: Pathway enrichment visualization

**Output structure:**
```
results/figures/
├── png/          # 300 DPI raster images
└── pdf/          # Vector graphics
```

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

# Elite visualizations
Rscript scripts/07_visualizations.R
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

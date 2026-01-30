# Bulk-RNA-Seq: A Step-by-Step Workflow (DESeq2 + QC + xCell + Enrichment)

A **highly documented, modular bulk RNA-seq workflow** using the **GEO breast cancer dataset GSE192341**.  
This repository is designed as a **scientifically rigorous reference pipeline**: each step explains **why it exists**, **what to check**, and **what we observed** in this dataset.

> **Note:** This repository focuses on **reproducible analysis logic + diagnostics**. Publication-grade figures and deep biological interpretation will be added in a future version.

---

## Purpose of This Repository

This repository provides a clear and reproducible pipeline for bulk RNA-seq differential expression analysis in R, covering:

1. **Data preparation from GEO processed matrices** (counts + gene annotation)
2. **Sample metadata retrieval from GEO** (building a DESeq2-ready `colData`)
3. **Input validation and gene filtering** (prevent silent mismatches)
4. **Differential expression analysis** with **DESeq2 + apeglm shrinkage**
5. **Dataset QC on transformed data** (PCA, sample distances, outlier flags)
6. **Cell-type deconvolution** with **xCell** (immune/stromal composition)
7. **Pathway enrichment** (Hallmark ORA/GSEA + GO BP + Reactome)

---

## Dataset

- **GEO Series:** [GSE192341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192341) (human breast tumors; neoadjuvant chemotherapy response)
- **Primary contrast:** **pCR (pathological complete response) vs No pCR**
- **Input:** GEO processed count-like matrix
  - `GSE192341_processed_data.txt.gz` (download from GEO supplementary files)

---

## Repository Structure

```
Bulk-RNA-Seq/
├── README.md
├── LICENSE
├── .gitignore
├── scripts/                          # Step-by-step pipeline scripts
│   ├── 00A_make_counts_from_processed_matrix.R
│   ├── 00B_make_coldata_from_GEO.R
│   ├── 01_input_validation.R
│   ├── 02_deseq2_analysis.R
│   ├── 03_qc_vst_pca_distances.R
│   ├── 04_deconvolution_xcell.R
│   └── 05_enrichment_ORA_GSEA.R
├── data/                             # Small files OK to commit
│   └── coldata.tsv
├── data_raw/                         # Large GEO downloads (DO NOT commit)
│   └── GSE192341_processed_data.txt
└── results/                          # Pipeline outputs (DO NOT commit)
    ├── tables/
    ├── figures/
    └── rds/
```

### What is tracked vs. not tracked
- ✅ **Tracked:** `scripts/`, `README.md`, `data/coldata.tsv`
- ❌ **NOT tracked:** `data_raw/`, `results/`, large matrices

---

## Installation & Setup

### Bioconductor Packages
```r
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
```

### CRAN Packages
```r
install.packages(c(
    "tidyverse",
    "data.table",
    "pheatmap",
    "RColorBrewer",
    "ggrepel",
    "msigdbr"
))

# xCell (from GitHub)
devtools::install_github("dviraran/xCell")
```

> **Tip:** In a future version, we will add `renv` for one-click reproducibility.

---

## How to Run (from repository root)

```bash
# Step 00A — Build counts matrix from GEO processed file
Rscript scripts/00A_make_counts_from_processed_matrix.R

# Step 00B — Build sample metadata (colData) from GEO
Rscript scripts/00B_make_coldata_from_GEO.R

# Step 01 — Input validation + QC + gene filtering
Rscript scripts/01_input_validation.R

# Step 02 — Differential expression (DESeq2 + apeglm shrinkage)
Rscript scripts/02_deseq2_analysis.R

# Step 03 — QC on VST-transformed data (PCA, distances, outliers)
Rscript scripts/03_qc_vst_pca_distances.R

# Step 04 — Cell-type deconvolution (xCell)
Rscript scripts/04_deconvolution_xcell.R

# Step 05 — Pathway enrichment (ORA + GSEA)
Rscript scripts/05_enrichment_ORA_GSEA.R
```

---

## Step-by-Step Workflow

Each section below explains: (1) **Purpose** — why this step exists, (2) **What to look for** — critical quality checkpoints, and (3) **What we found** — observations from GSE192341.

---

## Step 00A — Parse GEO Processed Matrix → Counts

### Purpose
GEO processed files often contain mixed data: gene annotations interleaved with expression values. This step:
- Converts the GEO processed file into a clean **gene × sample count matrix**
- Extracts a separate **gene annotation table** (Ensembl ID ↔ gene symbol) for downstream ID mapping

### What to look for
- Are sample columns correctly detected (numeric-like IDs)?
- Are gene IDs unique and valid Ensembl IDs (ENSG format)?
- Do matrix dimensions match expectations (~20,000+ genes, ~80+ samples)?

### What we found (GSE192341)
- ✅ Count matrix generated successfully with sample IDs (e.g., `2288`, `2292`, ...)
- ✅ Gene annotation table saved for future ENSEMBL → SYMBOL → ENTREZID mapping

### Key outputs
- `data/counts.csv`
- `data/gene_annotation.csv`

---

## Step 00B — Construct colData from GEO Metadata

### Purpose
DESeq2 requires a **sample metadata table** (`colData`) where:
- Row names exactly match count matrix column names
- A `condition` factor defines the contrast (here: `pCR` vs `No_pCR`)

This step queries GEOquery to extract clinical variables and build a DESeq2-ready `colData`.

### What to look for
- Do sample IDs in `colData` **exactly match** count matrix columns?
- What is the distribution of experimental groups?
- Are there samples with missing (`NA`) condition values?

### What we found (GSE192341)
- ✅ 87 samples loaded from GEO metadata
- ✅ Condition distribution: **No pCR = 61**, **pCR = 24**, **NA = 2**
- ⚠️ Two samples with missing response status — these are excluded from DE analysis

### Key outputs
- `data/coldata.tsv`

---

## Step 01 — Input Validation + Gene Filtering

### Purpose
**Sample mismatch is the most common silent error in bulk RNA-seq analysis.** This step:
- Validates that count matrix columns exactly match `colData` row names
- Computes basic QC metrics (library size, detected genes per sample)
- Applies **low-count gene filtering** to remove uninformative features

### Gene filtering rationale
Following DESeq2 best practices, we retain genes with **≥10 counts in ≥3 samples**. This:
- Reduces the multiple testing burden
- Removes genes with insufficient information for reliable dispersion estimation
- Improves statistical power for true DE genes

### What to look for
- ✅ Sample alignment check passes (no mismatch warnings)
- Library sizes: no extreme outliers (possible failed libraries or contamination)
- Genes retained after filtering: typically ~15,000–18,000 protein-coding genes remain

### What we found (GSE192341)
- ✅ Sample matching passed: count columns match `colData` rows
- ✅ Gene filtering applied: reduced from ~60,000 → ~18,000 features
- QC metrics exported for downstream visualization

### Key outputs
- `results/tables/sample_qc_metrics.tsv`
- `results/tables/filtered_counts.tsv`
- `results/tables/input_summary.tsv`

---

## Step 02 — Differential Expression (DESeq2 + apeglm)

### Purpose
Perform **negative binomial-based differential expression testing** using DESeq2:
- Estimate size factors (library depth normalization)
- Estimate gene-wise dispersions (variance modeling)
- Fit GLM and perform Wald test for condition effect
- Apply **apeglm shrinkage** to obtain regularized log2 fold changes

### Why apeglm shrinkage?
Raw MLE (maximum likelihood estimate) log2 fold changes can be noisy for low-count genes. The apeglm method applies an adaptive, heavy-tailed prior that:
- Shrinks noisy estimates toward zero
- Preserves large, well-supported fold changes
- Produces more stable rankings for downstream GSEA

### What to look for
- Number of significant genes at various thresholds (padj < 0.05, < 0.01)
- Balance of up- vs. down-regulated genes
- Sanity check: top DE genes should include known biology (e.g., immune markers in responders)

### What we found (GSE192341)
- ✅ DESeq2 completed successfully
- Significant genes (padj < 0.05): see `DESeq2_run_summary.tsv`
- Full results with shrunken LFC exported

### Key outputs
- `results/rds/dds.rds` — DESeqDataSet object
- `results/tables/DESeq2_results_all.tsv`
- `results/tables/DESeq2_results_significant_padj0.05.tsv`
- `results/tables/DESeq2_results_LFCshrink_apeglm.tsv`
- `results/tables/normalized_counts.tsv`
- `results/tables/vst_matrix.tsv`

---

## Step 03 — QC on Transformed Data (PCA, Distances, Outliers)

### Purpose
Before interpreting DE results, assess whether the **global transcriptomic structure** supports the experimental design:
- **PCA** reveals major sources of variation — do samples cluster by condition or by confounders?
- **Sample distance heatmap** highlights outliers and batch effects
- **Outlier detection** flags samples with unusually high mean pairwise distance

### Variance-stabilizing transformation (VST)
For visualization and distance calculations, raw counts must be transformed to stabilize variance across the dynamic range. VST is preferred over log2(counts+1) because it properly accounts for the mean-variance relationship in RNA-seq data.

### What to look for
- What percentage of variance is explained by PC1 and PC2?
- Do pCR and No_pCR samples show any separation, or overlap completely?
- Are there outlier samples that cluster separately (possible technical failure)?

### What we found (GSE192341)
- ✅ PCA coordinates and variance explained exported
- ✅ Outlier flag table generated (review samples with Z-score > 3)
- ⚠️ Limited separation between pCR/No_pCR on PC1–PC2 (common in heterogeneous tumor samples)

### Key outputs
- `results/tables/qc_pca_coordinates.tsv`
- `results/tables/qc_pca_variance_explained.tsv`
- `results/tables/qc_outlier_flags.tsv`
- `results/tables/qc_top_variable_genes.tsv`

---

## Step 04 — Cell-Type Deconvolution (xCell)

### Purpose
Bulk RNA-seq measures a mixture of cell types. **Deconvolution algorithms** estimate the relative abundance of different cell populations from the bulk expression profile. xCell provides enrichment scores for 64 immune and stromal cell types.

### Why deconvolution matters
In cancer immunotherapy/chemotherapy studies, treatment response often correlates with:
- CD8+ T cell infiltration (cytotoxic immune response)
- M1 vs M2 macrophage polarization
- Regulatory T cell (Treg) abundance

Comparing these scores between pCR and No_pCR may reveal microenvironmental signatures of response.

### What to look for
- Successful ENSEMBL → SYMBOL mapping (should convert thousands of genes)
- Cell score matrix dimensions: samples × 64 cell types
- Which cell types show significant differences (Wilcoxon test, BH-adjusted p-value)?

### What we found (GSE192341)
- ✅ ENSEMBL → SYMBOL mapping successful
- ✅ xCell score matrix generated
- Group comparison table exported for visualization (boxplots, heatmaps pending)

### Key outputs
- `results/tables/deconvolution_xcell_scores.tsv`
- `results/tables/deconvolution_xcell_pcr_vs_nopcr.tsv`
- `results/tables/deconvolution_xcell_run_summary.tsv`

---

## Step 05 — Pathway Enrichment (ORA + GSEA)

### Purpose
Individual gene-level results are difficult to interpret biologically. **Pathway enrichment** aggregates DE signals into functionally coherent gene sets:
- **ORA (Over-Representation Analysis):** Tests whether DE genes are enriched in a pathway beyond random chance
- **GSEA (Gene Set Enrichment Analysis):** Uses the full ranked gene list (no arbitrary cutoff) to detect coordinated shifts

### Gene sets used
| Collection | Source | Description |
|------------|--------|-------------|
| Hallmark | MSigDB | 50 well-defined biological processes |
| GO BP | Gene Ontology | Biological Process terms |
| Reactome | Reactome | Curated pathway database |

### What to look for
- Successful ENSEMBL → ENTREZID mapping (required for GO/Reactome)
- Do enriched terms make biological sense (immune response, cell cycle, EMT)?
- Consistency between ORA and GSEA results

### What we found (GSE192341)
- ✅ Annotated DE table generated (ENSEMBL + SYMBOL + ENTREZID)
- ✅ Enrichment outputs generated for Hallmark, GO BP, and Reactome
- Top pathways to be interpreted in future "Biological Interpretation" section

### Key outputs
- `results/tables/DESeq2_results_annotated.tsv`
- `results/tables/enrich_hallmark_ORA_UP.tsv`
- `results/tables/enrich_hallmark_ORA_DOWN.tsv`
- `results/tables/enrich_GO_BP_ORA_UP.tsv`
- `results/tables/enrich_GO_BP_ORA_DOWN.tsv`
- `results/tables/enrich_reactome_ORA_UP.tsv`
- `results/tables/enrich_reactome_ORA_DOWN.tsv`
- `results/tables/enrich_hallmark_GSEA.tsv`

---

## Limitations & Planned Improvements

- [ ] Add publication-grade visualizations (volcano plots, MA plots, PCA biplots, enrichment dotplots, xCell heatmaps)
- [ ] Incorporate covariates (molecular subtype, ER/HER2 status) into multifactor DESeq2 design
- [ ] Add batch effect diagnostics (if batch metadata becomes available)
- [ ] Generate Quarto/Rmd report with GitHub Pages deployment
- [ ] Add `renv` lockfile for full computational reproducibility

---

## References & Citations

If you use or adapt this workflow, please cite the original methods:

| Tool | Citation |
|------|----------|
| **DESeq2** | Love MI, Huber W, Anders S. *Genome Biology* 2014;15:550 |
| **apeglm** | Zhu A, Ibrahim JG, Love MI. *Bioinformatics* 2019;35:2084–2092 |
| **xCell** | Aran D, Hu Z, Butte AJ. *Genome Biology* 2017;18:220 |
| **clusterProfiler** | Wu T et al. *The Innovation* 2021;2:100141 |
| **MSigDB/Hallmark** | Liberzon A et al. *Cell Systems* 2015;1:417–425 |
| **GEO Dataset** | GSE192341 |

---

## Contact

📬 For questions, suggestions, or issues:

**Email:** bioinfosourabh@gmail.com  
**Website:** [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app)

---

## License

This project is open source and available under the MIT License.

# Bulk-RNA-Seq: A Step-by-Step Workflow (DESeq2 + QC + xCell + Enrichment)

A **highly documented, modular bulk RNA-seq workflow** using **GEO breast cancer dataset GSE192341**.  
This repository is designed as a **scientifically rigorous reference pipeline**: each step explains **why it exists**, **what to check**, and **what we observed** in this dataset.

> Note: This repo focuses on **reproducible analysis logic + diagnostics**. Publication-grade figures and deep biological interpretation will be added in a later version.

---

## Purpose of This Repository

This repository provides a clear and reproducible pipeline for bulk RNA-seq analysis in R, covering:

1. **Data preparation from GEO processed matrices** (counts + gene annotation)
2. **Sample metadata retrieval from GEO** (building a DESeq2-ready `colData`)
3. **Input validation and gene filtering** (prevent silent mismatches)
4. **Differential expression analysis** with **DESeq2 + apeglm shrinkage**
5. **Dataset QC on transformed data** (PCA, distances, outlier flags)
6. **Cell-type deconvolution** (xCell; immune/stromal composition)
7. **Pathway enrichment** (Hallmark ORA/GSEA + GO BP + Reactome)

---

## Dataset

- **GEO Series:** GSE192341 (human breast tumors; neoadjuvant response endpoints)
- **Primary contrast used here:** **pCR vs No pCR**
- Input is GEO *processed* count-like matrix:
  - `GSE192341_processed_data.txt.gz` (download from GEO supplementary files)

---

## Repository structure


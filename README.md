# Breast Cancer Gene Expression Analysis using TCGA-BRCA (RNA-seq workflow in R)

## Problem Statement:
   Gene expression analysis is a core tool in modern cancer research. While simplified "Tumor vs Normal" comparisons are common entry points, real-world biological insights require interrogating **molecular subtypes** and **functional pathways**.

   I created this project to demonstrate a rigid, end-to-end RNA-seq analysis pipeline applied to the TCGA-BRCA cohort. This project moves beyond basic tutorials by addressing clinical heterogeneity (Subtypes) and functional context (GSEA).

## Objectives:
   1. **Build a reproducible pipeline**: From raw count download to functional enrichment.
   2. **Interrogate Molecular Subtypes**: specifically comparing **Basal-like** (aggressive) vs **Luminal A** (hormone receptor-positive) tumors.
   3. **Identify Drivers**: Detect differential expressed genes defining these subtypes.
   4. **Pathway Analysis**: Use Gene Set Enrichment Analysis (GSEA) to understand the biological processes varying between subtypes.

## Data Source:
   The analysis uses the complete available set of Primary Tumor samples from the **TCGA Breast Invasive Carcinoma (TCGA-BRCA)** project.
   - **Data Fetching**: Automated via `TCGAbiolinks` to ensure reproducibility.
   - **Subtyping**: Samples are filtered using **PAM50** labels (Luminal A vs Basal).

## Methods Overview:
 **Data Preparation**:
   Downloads the full TCGA-BRCA RNA-seq cohort. Clinical metadata is merged with expression data to assign PAM50 subtypes.

 **Quality Control**:
   Library size assessment and Variance Stabilizing Transformation (VST) specifically for visualizing subtype separation.

 **Differential Expression (Basal vs Luminal A)**:
   Uses `DESeq2` to model count data: `~ BRCA_Subtype_PAM50`.
   - **Contrast**: Basal vs Luminal A.
   - **Statistical Rigor**: Corrects for multiple testing (Benjamini-Hochberg).

 **Gene Set Enrichment Analysis (GSEA)**:
   Performs functional enrichment using `fgsea` and the **MSigDB Hallmark** gene sets. This identifies upregulated pathways (e.g., Cell Cycle in Basal) vs downregulated ones (e.g., Estrogen Response).

 **Visualization**:
   - **PCA**: To visualize global separation of subtypes.
   - **Volcano Plots**: To view effect sizes vs significant.
   - **Heatmaps**: To view top gene signatures.
   - **GSEA Plots**: To view enriched pathways.

## Limitations & Future Work:
   1. **Survival Integration**: Future phases will integrate clinical survival data (Unix & Cox regression) to link gene expression to patient outcomes.
   2. **Multi-Omics**: Integrating CNV or Methylation data would provide a systems biology view.

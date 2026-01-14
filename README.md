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
  
## Quick Start (End-to-End)
1. **Install dependencies** (Bioconductor + CRAN):
   ```r
   source("scripts/install_dependencies.R")
   ```
2. **Run the full pipeline** (download → QC → DE → GSEA):
   ```r
   source("run_pipeline.R")
   ```

> **Note**: Data downloads can take time depending on network speed.

## Outputs (Generated Files)
**Results tables**
- `results/deseq2_results_Basal_vs_LuminalA.csv`: DESeq2 results.
- `results/gsea_results.csv`: GSEA results.
- `results/subtype_counts_raw.txt`: subtype counts before downsampling.
- `results/subtype_counts_downsampled.txt`: subtype counts after downsampling.
- `results/sessionInfo.txt`: captured R session info for reproducibility.

**Figures**
- `figures/pca_subtypes.png`
- `figures/volcano_basal_vs_luminalA.png`
- `figures/heatmap_top_genes.png`
- `figures/gsea_summary.png`
- `figures/qc_library_size.png`

## Reproducibility Notes
- **Fixed random seed**: used during downsampling to keep results consistent across runs.
- **Session tracking**: `results/sessionInfo.txt` is written at the end of the pipeline.
- **Deterministic inputs**: data is sourced programmatically from TCGA via `TCGAbiolinks`.

## Data & Design Choices
- **Cohort**: Primary Tumor samples from TCGA-BRCA.
- **Subtypes**: PAM50 labels filtered to **Basal** and **Luminal A**.
- **DE model**: `~ group` with Luminal A as the reference level.
- **Filtering**: genes with total counts < 10 are removed prior to DESeq2.

## Interpretation Pointers (What to Look For)
- **PCA**: subtype separation along PC1/PC2 indicates global expression differences.
- **Volcano plot**: highlights statistically significant genes with large effects.
- **GSEA**: pathways with high |NES| indicate subtype-specific biology.

## Limitations & Future Work:
   1. **Survival Integration**: Future phases will integrate clinical survival data (Unix & Cox regression) to link gene expression to patient outcomes.
   2. **Multi-Omics**: Integrating CNV or Methylation data would provide a systems biology view.
   3. **Covariates & batch**: Adding clinical covariates (age, stage, purity) or batch variables would reduce confounding.
   

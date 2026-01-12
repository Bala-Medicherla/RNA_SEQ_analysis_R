# Breast Cancer Gene Expression Analysis using TCGA-BRCA (RNA-seq workflow in R)

## Problem Statement:
   Gene expression analysis is a core tool in modern cancer research, but working with RNA-seq data involves much more than running a few       statistical tests. Before any biological conclusions can be made, analysts must carefully organize data, perform quality checks, choose      appropriate transformations, and explore patterns at both the sample and gene level.
   
   I created this project to practice and document the end to end thought process behind a basic RNA-seq analysis, using breast cancer as a     motivating example. The focus is on methodology, clarity, and reproducibility, rather than on discovering new biological findings.

## Objectives:
   The goal of this analysis are intentionally modest and educational:
   1. Build a clean, reproducible RNA-seq–style analysis workflow in R.
   2. Explore whether tumor and normal breast tissue samples show global expression differences.
   3. Identify genes that appear differentially expressed between groups in a simplified setting.
   4. Practice standard visualization techniques commonly used in transcriptomics.
   5. Clearly document assumptions, limitations, and next steps

## Data Source:
   The analysis uses publicly available breast cancer gene expression data derived from the TCGA Breast Invasive Carcinoma (TCGA-BRCA)          project.
   To keep the workflow lightweight and reproducible.
   1. A curated subset of tumor and normal samples is used.
   2. Only expression data and minimal sample annotations are included.
   3. No controlled or patient identifying information is accessed.
    
## Methods Overview:
 **Data preparation**:
   Expression values and sample metadata are loaded and aligned. Sample labels are checked to ensure that tumor and normal groups are           correctly defined before any analysis begins.
   
 **Quality Control**:
   Simple quality control checks are performed by examining library sizes (total expression per sample). These plots help confirm that no       samples show extreme behavior that could distort downstream analyses.
   
 **Transformation for exploration**:
   A log based transformation is applied to stabilize variance for visualization purposes. This step is used only for exploratory analysis      and is clearly distinguished from formal statistical modeling.
   
 **Differential expression(tumor vs normal)**:
   Genes are compared between tumor and normal samples to identify expression differences. Statistical testing is performed at the gene         level, and p-values are adjusted to control the false discovery rate. Results are interpreted cautiously, with attention to both effect      size and statistical significance.
   
 **Visualization**:
   Several standard plots are generated to summarize the data:
   1. Principal Component Analysis (PCA) to assess sample-level separation.
   2. Volcano plots to visualize gene-level effects and significance.
   3. Heatmaps to examine expression patterns among top-ranked genes.


## What the results show:
   At a high level, the analysis illustrates patterns commonly seen in breast cancer expression studies:
   1. Tumor and normal samples tend to separate in low-dimensional space.
   2. A subset of genes shows consistent expression differences between groups.
   3. Top differentially expressed genes display structured patterns across samples.
  These findings are illustrative, not confirmatory, and are used to practice interpretation rather than draw biological conclusions.

## Limitations:
   This project intentionally keeps scope limited:
   1. Only a subset of samples is analyzed.
   2. Differential expression is performed using a simplified approach.
   3. No pathway or gene-set enrichment analysis is included.
   These choices were made to prioritize clarity and learning over completeness.

## Honest note:
   This project reflects my effort to move beyond running code towards **understanding the reasoning behind genomic data analysis**. It         represents a learning step.

# Breast Cancer Gene Expression Analysis (Tumor vs Normal)
simple R-based analysis modeled after RNA-seq workflows, using simulated breast cancer data.

## Overview
In this project, I used R to carry out an analysis inspired by standard RNA-seq workflows, reflecting approaches commonly applied in biomedical informatics and cancer research.

## Objective:
This project focuses on developing a practical understanding of gene expression analysis, with an emphasis on working in R, visualizing biological data, and performing basic statistical analyses.

## Dataset:
The dataset used in this analyses is a simulated breast cancer gene expression matrix with:
3,000 genes
20 samples
10 Tumor
10 Normal
log2 fold-change patterns

The dataset includes two files:
expression_matrix – gene IDs, gene symbols, and expression counts
sample_metadata – sample IDs, tumor/normal labels, and subtype
This structure resembles public breast cancer datasets (e.g., TCGA-BRCA).

## Analysis Workflow:
This project follows a simplified but realistic genomics analysis workflow, similar to the early stages of an RNA-seq study.

  1. **Data import and preparation**:
     For this project, I generated a simulated RNA-seq dataset directly in R to reflect a realistic study setup. The gene expression values and sample information were      then organized into a clean matrix format so the data could be used for downstream analysis
     Basic quality checks were performed by calculating library sizes (total expression per sample). These were visualized using bar plots to ensure that no samples showed     extreme or unexpected behavior.

  3. **Log transformation**:
      A log₂(x + 1) transformation was applied to stabilize variance across genes. This step helps make downstream visualizations, such as PCA and heatmaps, easier to interpret.

  4. **Differential expression analysis (Tumor vs. Normal)**:
      Differential expression was assessed by performing a simple t-test for each gene to compare tumor and normal samples. For each gene, log₂ fold change, p-values, and FDR-adjusted p-values (using the Benjamini–Hochberg method) were calculated, and genes were ranked by statistical significance.

  5. **Visualization and interpretation**:
      Plots which are generated to summarize the results:
      Volcano plot to highlight up- and down-regulated genes.
      PCA plot to visualize sample-level separation between tumor and normal groups.
      Heatmap showing the top 50 differentially expressed genes.

## Key Results:

The PCA shows that tumor and normal samples cluster separately based on gene expression, with additional variation likely reflecting breast cancer subtypes. The overall pattern matches what is typically observed in breast cancer studies.
![image alt](https://github.com/Bala-Medicherla/Breast-cancer-Gene-expression-analysis/blob/9837cead6d7ec579d7669ec7c8afd45fb221c1a5/pca_tumor_normal.png)

Bar plot shows the total gene expression counts for each sample. Tumor and normal samples display fairly consistent total expression levels within their respective groups, with no samples standing out as extreme outliers. This suggests that sequencing depth (or overall expression magnitude in this simulated dataset) is comparable across samples.
![image alt](https://github.com/Bala-Medicherla/Breast-cancer-Gene-expression-analysis/blob/26da6c9ce24c9b0960b71be5b52adf6a08694542/total_counts_barplot.png)

Volcano plot showing differentially expressed genes between tumor and normal samples. Significant genes meet both fold-change and adjusted p-value thresholds, with many showing higher expression in tumors.
while the heatmap reveals structured expression patterns across samples.
![image alt](https://github.com/Bala-Medicherla/Breast-cancer-Gene-expression-analysis/blob/d71458737edceb03aa74328526df76029db184ff/volcano_plot_tumor_vs_normal.png)

Quality control checks indicate consistent sample depth across all samples.
Overall, these results are in line with patterns commonly reported in breast cancer gene expression studies.

## Skills Demonstrated:
This project helped me develop and demonstrate:
R programming fundamentals
Data wrangling with base R and tidy syntax
ggplot2 visualization
Differential expression concepts
Quality control for genomic datasets
Reproducible workflow design

## Next Steps
As I continue learning, I plan to:
Analyze real breast cancer datasets (TCGA-BRCA)
Explore pathway enrichment analysis (clusterProfiler)
Perform survival analysis using clinical variables

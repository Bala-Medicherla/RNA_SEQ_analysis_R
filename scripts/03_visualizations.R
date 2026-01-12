############################################################
# 03_visualizations.R
#
# Goal:
#   Generate standard RNA-seq visualizations that summarize
#   expression patterns and differential expression results.
#
# Why this step matters:
#   Visual summaries (PCA, volcano plots, heatmaps) help
#   assess global structure, effect sizes, and consistency
#   of gene-level signals in a way tables alone cannot.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

# ----------------------------------------------------------
# Load objects generated in previous steps
# ----------------------------------------------------------

vsd <- readRDS("data/processed/vsd_object.rds")
dds <- readRDS("data/processed/dds_object.rds")
de_results <- read.csv("results/deseq2_results_tumor_vs_normal.csv")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# PCA: sample-level structure
# ----------------------------------------------------------

# PCA helps assess whether tumor and normal samples
# separate based on global expression patterns
pca_df <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA of TCGA-BRCA Samples",
    x = "PC1",
    y = "PC2"
  )

ggsave(
  filename = "figures/pca_tumor_vs_normal.png",
  plot     = p_pca,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# ----------------------------------------------------------
# Volcano plot: effect size vs significance
# ----------------------------------------------------------

# Volcano plots summarize differential expression by
# combining fold change and statistical significance
p_volcano <- ggplot(
  de_results,
  aes(x = log2FoldChange, y = -log10(padj))
) +
  geom_point(alpha = 0.4) +
  labs(
    title = "Volcano Plot: Tumor vs Normal",
    x = "log2 Fold Change",
    y = "-log10(FDR-adjusted p-value)"
  )

ggsave(
  filename = "figures/volcano_tumor_vs_normal.png",
  plot     = p_volcano,
  width    = 6,
  height   = 5,
  dpi      = 300
)

# ----------------------------------------------------------
# Heatmap: top differentially expressed genes
# ----------------------------------------------------------

# Select top genes by adjusted p-value
top_genes <- de_results |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  slice_head(n = 30) |>
  pull(gene_id)

# Extract expression values for selected genes
heatmap_matrix <- assay(vsd)[top_genes, ]

pheatmap(
  heatmap_matrix,
  show_rownames = FALSE,
  filename = "figures/heatmap_top_genes.png"
)

message("All visualizations generated successfully.")
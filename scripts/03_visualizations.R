############################################################
# 03_visualizations.R
#
# Goal:
#   Visualize Basal vs Luminal A differences.
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(DESeq2)
})

# ----------------------------------------------------------
# Load objects
# ----------------------------------------------------------

vsd <- readRDS("data/processed/vsd_object.rds")
# dds <- readRDS("data/processed/dds_object.rds") # Not strictly needed if we have vsd & results
de_results <- read.csv("results/deseq2_results_Basal_vs_LuminalA.csv")

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# PCA
# ----------------------------------------------------------

pca_df <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA: Basal vs Luminal A Subtypes",
    subtitle = "Separation of breast cancer subtypes",
    x = "PC1",
    y = "PC2",
    color = "Subtype"
  )

ggsave("figures/pca_subtypes.png", p_pca, width = 6, height = 5)

# ----------------------------------------------------------
# Volcano Plot
# ----------------------------------------------------------

p_volcano <- ggplot(de_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Basal vs Luminal A",
    x = "log2 Fold Change (Basal / Luminal A)",
    y = "-log10(FDR)"
  ) +
  scale_color_manual(values = c("gray", "red")) +
  theme(legend.position = "none")

ggsave("figures/volcano_basal_vs_luminalA.png", p_volcano, width = 6, height = 5)

# ----------------------------------------------------------
# Heatmap: Top 30 Genes
# ----------------------------------------------------------

top_genes <- de_results |>
  filter(!is.na(padj)) |>
  arrange(padj) |>
  slice_head(n = 30) |>
  pull(gene_id)

heatmap_matrix <- assay(vsd)[top_genes, ]

# Add annotation for columns
annotation_col <- data.frame(
  Subtype = colData(vsd)$group
)
rownames(annotation_col) <- colnames(vsd)

pheatmap(
  heatmap_matrix,
  show_rownames = TRUE, # Show gene names for top 30
  show_colnames = FALSE,
  annotation_col = annotation_col,
  main = "Top 30 Differentially Expressed Genes",
  filename = "figures/heatmap_top_genes.png",
  scale = "row" # Z-score scaling makes patterns clearer
)

message(" visualizations generated.")
############################################################
# 04_gsea_analysis.R
#
# Goal:
#   Perform Gene Set Enrichment Analysis (GSEA) to understand
#   functions enriched in Basal vs Luminal A subtypes.
#
# Method:
#   - Use 'fgsea' for fast GSEA.
#   - Use 'msigdbr' to fetch Hallmark gene sets.
############################################################

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------------------------------------
# Load Data
# ----------------------------------------------------------

de_results <- read.csv("results/deseq2_results_Basal_vs_LuminalA.csv")

# ----------------------------------------------------------
# Prepare Ranked List
# ----------------------------------------------------------

# GSEA requires a ranked list of genes.
# Using log2FoldChange as the ranking metric for simplicity
# (Positive = Upregulated in Basal, Negative = Upregulated in Luminal A)
ranked_genes <- de_results |>
  filter(!is.na(log2FoldChange) & !is.na(padj)) |>
  arrange(desc(log2FoldChange)) |>
  select(gene_id, log2FoldChange)

ranks <- setNames(ranked_genes$log2FoldChange, ranked_genes$gene_id)

# ----------------------------------------------------------
# Get Gene Sets (Hallmark)
# ----------------------------------------------------------

# Fetch Hallmark gene sets for Homo sapiens
# Note: 'category' is used for compatibility. 
# If a warning appears suggesting 'collection', it is safe to ignore 
# as long as the gene sets are retrieved.
m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(x = m_df$gene_symbol, f = m_df$gs_name)

# ----------------------------------------------------------
# Run GSEA
# ----------------------------------------------------------

message("Running GSEA...")
fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize=15, maxSize=500)

# Sort by NES (Normalized Enrichment Score)
fgseaRes <- fgseaRes |> arrange(desc(NES))

# ----------------------------------------------------------
# Visualize
# ----------------------------------------------------------

# Top pathways
topPathways <- bind_rows(
  head(fgseaRes, 10),
  tail(fgseaRes, 10)
)

p_gsea <- ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(
    title = "GSEA: Hallmark Pathways (Basal vs Luminal A)",
    x = "Pathway",
    y = "Normalized Enrichment Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("figures/gsea_summary.png", p_gsea, width = 8, height = 6)

# Save tables
write.csv(as.data.frame(fgseaRes), "results/gsea_results.csv", row.names = FALSE)

message("GSEA analysis completed.")

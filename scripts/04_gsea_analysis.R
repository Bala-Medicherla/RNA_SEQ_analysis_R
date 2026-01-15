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
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ----------------------------------------------------------
# Load Data
# ----------------------------------------------------------

de_results <- read.csv("results/deseq2_results_Basal_vs_LuminalA.csv")

# ----------------------------------------------------------
# Prepare Ranked List
# ----------------------------------------------------------

# GSEA requires a ranked list of genes (symbols for MSigDB).
# Using log2FoldChange as the ranking metric for simplicity
# (Positive = Upregulated in Basal, Negative = Upregulated in Luminal A)
ranked_genes <- de_results |>
  dplyr::filter(!is.na(log2FoldChange) & !is.na(padj)) |>
  dplyr::arrange(desc(log2FoldChange)) |>
  dplyr::select(gene_id, log2FoldChange)

gene_ids <- ranked_genes$gene_id
if (any(grepl("^ENSG", gene_ids))) {
  message("Detected Ensembl IDs. Mapping to gene symbols for GSEA.")
  gene_ids <- gsub("\\..*$", "", gene_ids)
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(gene_ids),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  )
  ranked_genes <- ranked_genes |>
    dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::left_join(mapping, by = c("gene_id" = "ENSEMBL")) |>
    dplyr::filter(!is.na(SYMBOL)) |>
    dplyr::mutate(gene_id = SYMBOL) |>
    dplyr::select(gene_id, log2FoldChange)
}

ranked_genes <- ranked_genes |>
  dplyr::group_by(gene_id) |>
  dplyr::summarize(log2FoldChange = max(log2FoldChange), .groups = "drop") |>
  dplyr::arrange(desc(log2FoldChange))

if (nrow(ranked_genes) == 0) {
  stop("No genes available for ranking after filtering/mapping.")
}

gene_ids <- ranked_genes$gene_id
if (any(grepl("^ENSG", gene_ids))) {
  message("Detected Ensembl IDs. Mapping to gene symbols for GSEA.")
  gene_ids <- gsub("\\..*$", "", gene_ids)
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(gene_ids),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  )
  ranked_genes <- ranked_genes |>
    mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    left_join(mapping, by = c("gene_id" = "ENSEMBL")) |>
    filter(!is.na(SYMBOL)) |>
    mutate(gene_id = SYMBOL) |>
    select(gene_id, log2FoldChange)
}

ranked_genes <- ranked_genes |>
  group_by(gene_id) |>
  summarize(log2FoldChange = max(log2FoldChange), .groups = "drop") |>
  arrange(desc(log2FoldChange))

if (nrow(ranked_genes) == 0) {
  stop("No genes available for ranking after filtering/mapping.")
}

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
fgseaRes <- fgseaRes |> dplyr::arrange(desc(NES))

# ----------------------------------------------------------
# Visualize
# ----------------------------------------------------------

# Top pathways
if (nrow(fgseaRes) == 0) {
  stop("fgsea returned no pathways. Check gene ID mapping and rank list size.")
}

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

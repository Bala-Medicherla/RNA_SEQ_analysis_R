############################################################
# 02_deseq2_analysis.R
#
# Goal:
#   Identify genes that are differentially expressed between
#   tumor and normal breast tissue using DESeq2.
#
# Why this step matters:
#   RNA-seq data are counts, not continuous measurements.
#   DESeq2 models these counts using a negative binomial
#   framework that properly accounts for biological and
#   technical variability.
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

# Create results directory if needed
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Load data prepared in previous steps
# ----------------------------------------------------------

counts_matrix   <- readRDS("data/processed/counts_matrix.rds")
sample_metadata <- readRDS("data/processed/sample_metadata.rds")

# Identify sample type column defensively
sample_type_col <- grep("sample_type", names(sample_metadata), value = TRUE)[1]

if (is.na(sample_type_col)) {
  stop("Sample type column not found in sample metadata.")
}

# Define analysis groups
sample_metadata$group <- factor(sample_metadata[[sample_type_col]])

# Keep only Tumor vs Normal samples
keep_samples <- sample_metadata$group %in%
  c("Primary Tumor", "Solid Tissue Normal")

counts_matrix   <- counts_matrix[, keep_samples]
sample_metadata <- sample_metadata[keep_samples, ]

# ----------------------------------------------------------
# Construct DESeq2 dataset
# ----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData   = sample_metadata,
  design    = ~ group
)

# Remove genes with extremely low counts across all samples
dds <- dds[rowSums(counts(dds)) >= 10, ]

# ----------------------------------------------------------
# Run DESeq2 model
# ----------------------------------------------------------

dds <- DESeq(dds)

# Extract results for Tumor vs Normal comparison
de_results <- results(
  dds,
  contrast = c("group", "Primary Tumor", "Solid Tissue Normal")
)

# Convert to tidy data frame for inspection and export
de_results_df <- as.data.frame(de_results) |>
  tibble::rownames_to_column("gene_id") |>
  arrange(padj)

# ----------------------------------------------------------
# Save outputs for reporting and visualization
# ----------------------------------------------------------

write.csv(
  de_results_df,
  "results/deseq2_results_tumor_vs_normal.csv",
  row.names = FALSE
)

# Save DESeq2 object for downstream visualization
saveRDS(dds, "data/processed/dds_object.rds")

message("DESeq2 differential expression analysis completed.")
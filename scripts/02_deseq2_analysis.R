############################################################
# 02_deseq2_analysis.R
#
# Goal:
#   Identify genes that are differentially expressed between
#   Basal-like and Luminal A breast cancer subtypes.
#
# Changes:
#   - Contrast: Basal vs Luminal A (Reference = Luminal A)
#   - Design: ~ BRCA_Subtype_PAM50
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

# Ensure subtype column exists
subtype_col <- "BRCA_Subtype_PAM50"

if (!subtype_col %in% names(sample_metadata)) {
  stop("Subtype column not found in sample metadata.")
}

# Define analysis groups
sample_metadata$group <- factor(sample_metadata[[subtype_col]])

# Explicitly set "Luminal A" as the reference level
# valid levels: "Luminal A", "Basal"
sample_metadata$group <- relevel(sample_metadata$group, ref = "Luminal A")

message(paste("Comparison levels:", paste(levels(sample_metadata$group), collapse=" vs ")))

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

message("Running DESeq2 (this may take a moment)...")
dds <- DESeq(dds)

# Extract results for Basal vs Luminal A
# Since we set Luminal A as reference, the coefficient is "group_Basal_vs_Reference"
de_results <- results(dds, contrast = c("group", "Basal", "Luminal A"))

# Convert to tidy data frame for inspection and export
de_results_df <- as.data.frame(de_results) |>
  tibble::rownames_to_column("gene_id") |>
  arrange(padj)

# ----------------------------------------------------------
# Save outputs for reporting and visualization
# ----------------------------------------------------------

write.csv(
  de_results_df,
  "results/deseq2_results_Basal_vs_LuminalA.csv",
  row.names = FALSE
)

# Save DESeq2 object for downstream visualization
saveRDS(dds, "data/processed/dds_object.rds")

message("DESeq2 analysis (Basal vs Luminal A) completed.")
############################################################
# 00_prepare_tcga_data.R
#
# Goal:
#   Download and prepare TCGA-BRCA RNA-seq data, specifically
#   targeting the "Luminal A" vs "Basal" biological comparison.
#
# Changes in Phase 2 (Optimized):
#   - Fetches full cohort metadata.
#   - Filters query for Luminal A / Basal *BEFORE* loading data.
#   - Caps at 200 samples/group to fit in 16GB RAM.
############################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
})

# ----------------------------------------------------------
# Project-level configuration
# ----------------------------------------------------------

project_id <- "TCGA-BRCA"
workflow_type <- "STAR - Counts"
sample_types <- c("Primary Tumor")
max_samples_per_type <- 200

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Query TCGA via GDC
# ----------------------------------------------------------

message("Querying GDC for TCGA-BRCA RNA-seq data (Primary Tumor)...")

query <- GDCquery(
  project       = project_id,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = workflow_type,
  sample.type   = sample_types
)

# ----------------------------------------------------------
# 2. Filter Query by Subtype (RAM Optimization)
# ----------------------------------------------------------

message("Retrieving PAM50 subtype data to filter BEFORE download/load...")
subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")

# Extract metadata contents from the query object
query_results <- query$results[[1]]

# The query results have 'cases' (submitter_id) or 'sample.submitter_id'
# We need to map these to the 'patient' column in subtypes.
# 'cases' looks like: TCGA-AR-A1AS-01A-11R-A128-07
# 'patient' looks like: TCGA-AR-A1AS

# Add patient ID to query results
query_results$patient <- substr(query_results$cases, 1, 12)

# Merge with subtype info
# using all.x=TRUE to keep query rows, but we only want matches
query_merged <- merge(query_results, subtypes, by = "patient", all.x = FALSE)

# Handle column name variations
if ("paper_BRCA_Subtype_PAM50" %in% names(query_merged)) {
  query_merged$BRCA_Subtype_PAM50 <- query_merged$paper_BRCA_Subtype_PAM50
}

target_subtypes <- c("Luminal A", "Basal")

# Filter for targets
query_filtered <- query_merged[query_merged$BRCA_Subtype_PAM50 %in% target_subtypes, ]

# ----------------------------------------------------------
# 3. Downsample to avoid Memory Crash
# ----------------------------------------------------------

message("Downsampling to avoid 16GB RAM limit...")
set.seed(42)

# Split by subtype
indices_by_type <- split(seq_len(nrow(query_filtered)), query_filtered$BRCA_Subtype_PAM50)

# Select random indices
keep_indices <- unlist(lapply(indices_by_type, function(idx) {
  if (length(idx) > max_samples_per_type) {
    sample(idx, max_samples_per_type)
  } else {
    idx
  }
}))

query_final_df <- query_filtered[keep_indices, ]

# IMPORTANT: sync the 'results' slot of the query object
# We must use the original columns matching what GDCquery produced
original_cols <- colnames(query$results[[1]])
# We need to ensure we only keep columns that were in the original query result
# (merge added many subtype columns)

# However, we need to pass the FILTERED rows back to the query object
# match barcodes (cases)
final_cases <- query_final_df$cases

# Filter the original query object
query$results[[1]] <- query$results[[1]][query$results[[1]]$cases %in% final_cases, ]

message(paste("Final query size:", nrow(query$results[[1]]), "samples."))

# ----------------------------------------------------------
# 4. Download and Prepare (Now Safe)
# ----------------------------------------------------------

message("Downloading filtered RNA-seq files...")
GDCdownload(query, directory = "data/raw")

message("Preparing SummarizedExperiment object (Memory Safe)...")
se <- GDCprepare(query, directory = "data/raw")

# ----------------------------------------------------------
# 5. Re-attach Metadata (Safety Step)
# ----------------------------------------------------------

# While we filtered by subtype, the 'se' object colData might be raw GDC metadata.
# We need to attach the PAM50 labels definitively for downstream scripts.

# Get the patient IDs again from the SE object
se_meta <- as.data.frame(colData(se))
se_meta$patient <- substr(se_meta$barcode, 1, 12)

# Merge again (safe)
se_meta_merged <- merge(se_meta, subtypes, by = "patient", all.x = TRUE)

if ("paper_BRCA_Subtype_PAM50" %in% names(se_meta_merged)) {
  se_meta_merged$BRCA_Subtype_PAM50 <- se_meta_merged$paper_BRCA_Subtype_PAM50
}

# Ensure row order
rownames(se_meta_merged) <- se_meta_merged$barcode
se_meta_merged <- se_meta_merged[colnames(se), ]

# Final check
counts_matrix <- assay(se)
sample_metadata <- se_meta_merged

# ----------------------------------------------------------
# 6. Save
# ----------------------------------------------------------

saveRDS(counts_matrix,   "data/processed/counts_matrix.rds")
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

message("Saved processed TCGA objects (Optimized).")
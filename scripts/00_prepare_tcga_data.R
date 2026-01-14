############################################################
# 00_prepare_tcga_data.R (Debug Version)
#
# Goal:
#   Download and prepare TCGA-BRCA RNA-seq data.
#   DEBUG: Ensure both Basal and Luminal A are retrieved.
############################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
})

# ----------------------------------------------------------
# Configuration
# ----------------------------------------------------------

project_id <- "TCGA-BRCA"
workflow_type <- "STAR - Counts"
sample_types <- c("Primary Tumor")
max_samples_per_type <- 200

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# 1. Query & Filter (Memory Safe)
# ----------------------------------------------------------

message("Querying GDC...")
query <- GDCquery(
  project       = project_id,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = workflow_type,
  sample.type   = sample_types
)

message("Retrieving Subtypes...")
subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")

query_results <- query$results[[1]]
query_results$patient <- substr(query_results$cases, 1, 12)

# Merge
query_merged <- merge(query_results, subtypes, by = "patient", all.x = FALSE)

# Merge
query_merged <- merge(query_results, subtypes, by = "patient", all.x = FALSE)

# Note: TCGAquery_subtype returns "LumA", "Basal", "Her2", "LumB", "Normal"
# We do not need to rename 'paper_' columns as the default is 'BRCA_Subtype_PAM50'

target_subtypes <- c("LumA", "Basal")
query_filtered <- query_merged[query_merged$BRCA_Subtype_PAM50 %in% target_subtypes, ]

# Normalize labels for downstream readability ("LumA" -> "Luminal A")
query_filtered$BRCA_Subtype_PAM50 <- ifelse(
  query_filtered$BRCA_Subtype_PAM50 == "LumA", "Luminal A", query_filtered$BRCA_Subtype_PAM50
)

# Update target to match the new readable labels
target_subtypes <- c("Luminal A", "Basal")


# DEBUG: Print found subtypes
print("Subtypes found in query:")
print(table(query_filtered$BRCA_Subtype_PAM50))

if (length(unique(query_filtered$BRCA_Subtype_PAM50)) < 2) {
  stop("CRITICAL ERROR: Found less than 2 subtypes! Cannot perform comparison.")
}

# ----------------------------------------------------------
# 2. Downsample
# ----------------------------------------------------------

message("Downsampling...")
set.seed(42)
indices_by_type <- split(seq_len(nrow(query_filtered)), query_filtered$BRCA_Subtype_PAM50)

keep_indices <- unlist(lapply(indices_by_type, function(idx) {
  if (length(idx) > max_samples_per_type) {
    sample(idx, max_samples_per_type)
  } else {
    idx
  }
}))

query_final_df <- query_filtered[keep_indices, ]
final_cases <- query_final_df$cases

# Update Query Object
query$results[[1]] <- query$results[[1]][query$results[[1]]$cases %in% final_cases, ]

message(paste("Final query size:", nrow(query$results[[1]]), "samples."))

# ----------------------------------------------------------
# 3. Download & Prepare
# ----------------------------------------------------------

message("Downloading (Filtered)...")
GDCdownload(query, directory = "data/raw")

message("Preparing SE...")
se <- GDCprepare(query, directory = "data/raw")

# ----------------------------------------------------------
# 4. Attach Metadata (Double Check)
# ----------------------------------------------------------

se_meta <- as.data.frame(colData(se))
se_meta$patient <- substr(se_meta$barcode, 1, 12)

# Re-merge to ensure metadata is attached to the final SE object
se_meta_merged <- merge(se_meta, subtypes, by = "patient", all.x = TRUE)

# Normalize "LumA" -> "Luminal A" in the final object too
se_meta_merged$BRCA_Subtype_PAM50 <- ifelse(
  se_meta_merged$BRCA_Subtype_PAM50 == "LumA", "Luminal A", se_meta_merged$BRCA_Subtype_PAM50
)

# Ensure checking logic works
# Some samples might be NA if they were in query but merge failed?
se_meta_merged <- se_meta_merged[!is.na(se_meta_merged$BRCA_Subtype_PAM50), ]

# Restore order
rownames(se_meta_merged) <- se_meta_merged$barcode
se_meta_merged <- se_meta_merged[colnames(se), ] # This might introduce NAs if barcode missing? 
# Better: subset SE to match meta
common_barcodes <- intersect(colnames(se), rownames(se_meta_merged))
se <- se[, common_barcodes]
se_meta_merged <- se_meta_merged[common_barcodes, ]

counts_matrix <- assay(se)
sample_metadata <- se_meta_merged

# DEBUG: Print final table
print("Final Metadata Subtypes:")
print(table(sample_metadata$BRCA_Subtype_PAM50))

# ----------------------------------------------------------
# 5. Save
# ----------------------------------------------------------

saveRDS(counts_matrix,   "data/processed/counts_matrix.rds")
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

message("Saved.")
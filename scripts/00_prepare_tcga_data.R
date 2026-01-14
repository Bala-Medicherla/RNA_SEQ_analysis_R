############################################################
# 00_prepare_tcga_data.R
#
# Goal:
#   Download and prepare TCGA-BRCA RNA-seq data, specifically
#   targeting the "Luminal A" vs "Basal" biological comparison.
#
# Changes in Phase 2:
#   - REMOVED sample limit (analyzing full relevant cohort)
#   - ADDED PAM50 subtype filtering
#   - Merged clinical subtype labels with expression data
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

# We focus on Primary Tumor samples and then filter by Subtype
sample_types <- c("Primary Tumor")

# Output directories
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
# 2. Download and Prepare Expression Data
# ----------------------------------------------------------

message("Downloading RNA-seq files (this may take significant time for the full cohort)...")
# We will use 'method = "api"' which is standard.
GDCdownload(query, directory = "data/raw")

message("Preparing SummarizedExperiment object...")
se <- GDCprepare(query, directory = "data/raw")

# ----------------------------------------------------------
# 3. Retrieve and Merge Subtype Information (PAM50)
# ----------------------------------------------------------

message("Retrieving PAM50 subtype data...")
# Recent TCGAbiolinks versions use PanCancerAtlas_subtypes or similar
subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")

# Extract barcode structure to merge
# se colData has 'patient' column usually, or we use the barcode.
# The 'patient' column in 'se' matches 'patient' in 'subtypes'.

meta <- as.data.frame(colData(se))

# Merge metadata with subtypes
# "patient" is the common key (e.g., TCGA-A2-A04N)
meta_with_subtypes <- merge(meta, subtypes, by = "patient", all.x = TRUE)

# Restore row order to match the counts matrix columns
rownames(meta_with_subtypes) <- meta_with_subtypes$barcode
meta_with_subtypes <- meta_with_subtypes[colnames(se), ]

# ----------------------------------------------------------
# 4. Filter for Luminal A and Basal
# ----------------------------------------------------------

target_subtypes <- c("Luminal A", "Basal")

# The column is usually 'BRCA_Subtype_PAM50'
if (!"BRCA_Subtype_PAM50" %in% names(meta_with_subtypes)) {
  warning("PAM50 column not found. Checking available columns...")
  # Fallback logic could go here, but for now we assume standard structure
}

keep_mask <- meta_with_subtypes$BRCA_Subtype_PAM50 %in% target_subtypes

# Apply filter
counts_matrix <- assay(se)[, keep_mask]
sample_metadata <- meta_with_subtypes[keep_mask, ]

message(paste("Filtered to", nrow(sample_metadata), "samples."))
message(paste("Breakdown:", paste(table(sample_metadata$BRCA_Subtype_PAM50), collapse=" / ")))

# ----------------------------------------------------------
# 5. Save Clean Objects
# ----------------------------------------------------------

saveRDS(counts_matrix,   "data/processed/counts_matrix.rds")
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

message("Saved processed TCGA objects (Luminal A vs Basal).")
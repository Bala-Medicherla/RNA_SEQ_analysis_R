############################################################
# 00_prepare_tcga_data.R
#
# Goal:
#   Download and prepare a reproducible subset of real-world
#   RNA-seq gene expression data from the TCGA-BRCA project.
#
# Why this script exists:
#   In real genomics projects, data acquisition is separated
#   from analysis. This script focuses ONLY on obtaining
#   expression data and saving clean objects that downstream
#   scripts can reuse without re-downloading TCGA data.
#
# What this script does:
#   1. Queries the Genomic Data Commons (GDC) for TCGA-BRCA
#      RNA-seq STAR-count data
#   2. Selects Primary Tumor and Solid Tissue Normal samples
#   3. Limits sample size to keep the project lightweight
#   4. Saves gene-level counts and sample metadata as .rds
############################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# ----------------------------------------------------------
# Project-level configuration
# ----------------------------------------------------------

# TCGA project identifier
project_id <- "TCGA-BRCA"

# RNA-seq processing workflow used by GDC
workflow_type <- "STAR - Counts"

# Sample types used for a simple, interpretable comparison
sample_types <- c("Primary Tumor", "Solid Tissue Normal")

# Limit samples per group to keep runtime reasonable
max_samples_per_group <- 25

# Output directories
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Query TCGA via the Genomic Data Commons
# ----------------------------------------------------------

message("Querying GDC for TCGA-BRCA RNA-seq data...")

query <- GDCquery(
  project       = project_id,
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = workflow_type,
  sample.type   = sample_types
)

# ----------------------------------------------------------
# Select a reproducible subset of samples
# ----------------------------------------------------------

# TCGA metadata schemas can vary slightly over time,
# so we defensively search for a column describing sample type
query_metadata <- query$results[[1]]
sample_type_column <- grep("sample_type", names(query_metadata), value = TRUE)[1]

if (is.na(sample_type_column)) {
  stop("Could not identify sample type column in TCGA metadata.")
}

# Select up to N samples per group to keep analysis manageable
keep_indices <- unlist(lapply(sample_types, function(type) {
  indices <- which(query_metadata[[sample_type_column]] == type)
  indices[seq_len(min(length(indices), max_samples_per_group))]
}))

query$results[[1]] <- query_metadata[keep_indices, , drop = FALSE]

message("Selected a limited, reproducible subset of samples.")

# ----------------------------------------------------------
# Download and prepare expression data
# ----------------------------------------------------------

message("Downloading RNA-seq files from GDC (this may take time)...")
GDCdownload(query, directory = "data/raw")

message("Preparing expression data into a SummarizedExperiment object...")
se <- GDCprepare(query, directory = "data/raw")

# ----------------------------------------------------------
# Extract counts and metadata for downstream analysis
# ----------------------------------------------------------

# Gene-level raw counts (rows = genes, columns = samples)
counts_matrix <- assay(se)

# Sample-level metadata (phenotype, sample type, etc.)
sample_metadata <- as.data.frame(colData(se))

# ----------------------------------------------------------
# Save clean objects for downstream scripts
# ----------------------------------------------------------

saveRDS(counts_matrix,   "data/processed/counts_matrix.rds")
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

message("Saved processed TCGA objects:")
message("- data/processed/counts_matrix.rds")
message("- data/processed/sample_metadata.rds")

message("TCGA-BRCA data preparation step completed successfully.")
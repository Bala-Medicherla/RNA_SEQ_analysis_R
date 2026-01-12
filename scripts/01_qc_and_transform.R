############################################################
# 01_qc_and_transform.R
#
# Goal:
#   Perform basic quality control (QC) and transformation
#   of TCGA-BRCA RNA-seq data before differential expression.
#
# Why this step matters:
#   QC helps ensure that downstream results are not driven
#   by extreme samples or technical artifacts. Transformation
#   (VST) is used ONLY for visualization, not for modeling.
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Create output directory for figures if it does not exist
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------
# Load processed data from Script 00
# ----------------------------------------------------------

counts_matrix   <- readRDS("data/processed/counts_matrix.rds")
sample_metadata <- readRDS("data/processed/sample_metadata.rds")

# Identify the column describing sample type (defensive coding)
sample_type_col <- grep("sample_type", names(sample_metadata), value = TRUE)[1]

if (is.na(sample_type_col)) {
  stop("Sample type column not found in sample metadata.")
}

# Define analysis groups
sample_metadata$group <- factor(sample_metadata[[sample_type_col]])

# ----------------------------------------------------------
# Construct DESeq2 object (needed even for QC)
# ----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData   = sample_metadata,
  design    = ~ group
)

# ----------------------------------------------------------
# QC: Library size per sample
# ----------------------------------------------------------

# Total counts per sample give a quick check for
# extreme sequencing depth differences
library_sizes <- colSums(counts(dds))

qc_df <- data.frame(
  sample = names(library_sizes),
  library_size = as.numeric(library_sizes),
  group = sample_metadata$group
)

qc_plot <- ggplot(
  qc_df,
  aes(x = reorder(sample, library_size),
      y = library_size,
      fill = group)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "QC: Library Size per Sample",
    x = "Sample",
    y = "Total Read Counts"
  )

ggsave(
  filename = "figures/qc_library_size.png",
  plot     = qc_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# ----------------------------------------------------------
# Variance-stabilizing transformation (for visualization)
# ----------------------------------------------------------

# Size factors are estimated to normalize sequencing depth
dds <- estimateSizeFactors(dds)

# VST makes expression values more comparable across genes
# and samples for PCA and heatmaps
vsd <- vst(dds, blind = TRUE)

# Save transformed object for downstream scripts
saveRDS(vsd, "data/processed/vsd_object.rds")

message("QC and transformation step completed successfully.")
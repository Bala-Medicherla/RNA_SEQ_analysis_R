############################################################
# install_dependencies.R
#
# Goal: Run this script ONCE to install all required packages
#       for the RNA-seq pipeline.
############################################################

# 1. Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. List of required packages (CRAN + Bioconductor)
required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "dplyr",
  "ggplot2",
  "pheatmap",
  "fgsea",    # For GSEA
  "msigdbr",  # For Gene Sets
  "AnnotationDbi",
  "org.Hs.eg.db"
)

# 3. Install missing packages
message("Checking for installed packages...")
installed_pkgs <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed_pkgs) {
    message(paste("Installing:", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message(paste("Already installed:", pkg))
  }
}

message("All dependencies installed!")

message("Running TCGA-BRCA RNA-seq pipeline...")

source("scripts/00_prepare_tcga_data.R")
source("scripts/01_qc_and_transform.R")
source("scripts/02_deseq2_analysis.R")
source("scripts/03_visualizations.R")
source("scripts/install_dependencies.R")
source("scripts/04_gsea_analysis.R")


message("Pipeline completed. Check the figures/ and results/ folders.")
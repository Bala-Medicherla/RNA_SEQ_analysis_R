
#rm(list=ls())
install.packages(c( "dplyr", "ggplot2", "pheatmap"))


# Breast Cancer Expression Analysis (Tumor vs Normal)
# Load libraries
library(dplyr)  
library(ggplot2) 
library(pheatmap)


set.seed(123)          
n_genes <- 6000       
n_normal <- 12        
n_tumor  <- 12         
subtypes <- c("LumA", "LumB", "HER2", "Basal") 
prop_de  <- 0.08      
lfc_mean <- 1.3   
lfc_sd   <- 0.4  


# Sample IDs
sample_ids <- sprintf("SAMPLE_%02d", 1:(n_normal + n_tumor))

# Metadata
meta_df <- data.frame(
  SampleID   = sample_ids,
  Condition  = c(rep("Normal", n_normal), rep("Tumor", n_tumor)),
  Subtype    = c(rep("NormalLike", n_normal),
                 sample(subtypes, n_tumor, replace = TRUE)),
  stringsAsFactors = FALSE
)

# Gene info
gene_ids <- sprintf("GENE_%05d", 1:n_genes)
gene_symbols <- paste0("G", 1:n_genes)

gene_info <- data.frame(
  GeneID = gene_ids,
  GeneSymbol = gene_symbols,
  stringsAsFactors = FALSE
)

# Simulate expression counts
# simulating RNA-seq-like counts with negative binomial
# Base mean per gene (some genes low, some high)
base_mean <- 2^(rnorm(n_genes, mean = 6, sd = 1.2)) 
dispersion <- 0.25                                   

# Choose DE genes
n_de <- round(n_genes * prop_de)
de_idx <- sample(seq_len(n_genes), n_de)

# Assigning tumor log2FC for DE genes (half up, half down)
tumor_log2fc <- rnorm(n_de, mean = lfc_mean, sd = lfc_sd)
sign_flip <- sample(c(-1, 1), n_de, replace = TRUE)
tumor_log2fc <- tumor_log2fc * sign_flip

# Building mean matrix (genes x samples)
mu <- matrix(base_mean, nrow = n_genes, ncol = n_normal + n_tumor)
colnames(mu) <- sample_ids
rownames(mu) <- gene_ids

# Applying tumor effects to DE genes (only tumor columns)
tumor_cols <- (n_normal + 1):(n_normal + n_tumor)
mu[de_idx, tumor_cols] <- mu[de_idx, tumor_cols] * (2^(tumor_log2fc))

# Sample-specific library size factors (to mimic sequencing depth variability)
lib_factors <- runif(n_normal + n_tumor, 0.8, 1.25) 
mu <- sweep(mu, 2, lib_factors, "*")

# Generating counts from NB:
# In rnbinom: size = 1/dispersion is a common parameterization
size <- 1 / dispersion
expr_counts <- matrix(
  rnbinom(n_genes * (n_normal + n_tumor), mu = as.vector(mu), size = size),
  nrow = n_genes
)

rownames(expr_counts) <- gene_ids
colnames(expr_counts) <- sample_ids

# Building expr_df to match original structure (GeneID, GeneSymbol, then samples)
expr_df <- data.frame(
  GeneID = gene_info$GeneID,
  GeneSymbol = gene_info$GeneSymbol,
  expr_counts,
  check.names = FALSE
)

#  Preparing expression matrix

gene_info <- expr_df[, c("GeneID", "GeneSymbol")]
expr_mat  <- expr_df[, -(1:2)]  # drop GeneID and GeneSymbol

# Converting to numeric matrix
expr_mat <- as.matrix(expr_mat)
str(expr_mat) 
is.numeric(expr_mat[, "SAMPLE_01"])
apply(expr_mat , 2,is.numeric) 


# Making sure column names match SampleID in metadata
colnames(expr_mat) <- colnames(expr_df)[-(1:2)]
stopifnot(all(colnames(expr_mat) == meta_df$SampleID))


# Basic QC: library sizes
# Total expression (like library size) per sample

total_counts <- colSums(expr_mat)

qc_df <- data.frame(
  SampleID = colnames(expr_mat),
  TotalCounts = total_counts,
  Condition = meta_df$Condition
)

# Barplot of total counts
p_qc <- ggplot(qc_df, aes(x = SampleID, y = TotalCounts, fill = Condition)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 2)) +
  ggtitle("Total Expression per Sample")

print(p_qc)

# Save QC plot
if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/total_counts_barplot.png", p_qc, width = 9, height = 4, dpi = 300)


#  Log2 transform
# Adding 1 to avoid log(0) , using log to become data easier(to reduce extreme values) to visualize and analyze
log_expr <- log2(expr_mat + 1) 

# Tumor vs Normal differential expression (simple t-test)

condition <- factor(meta_df$Condition, levels = c("Normal", "Tumor"))

# Using this Function to run t-test for each gene
t_test_results <- apply(log_expr, 1, function(x) {
  t.test(x[condition == "Tumor"],
         x[condition == "Normal"])$p.value
})

# Mean expression per group
mean_tumor  <- rowMeans(log_expr[, condition == "Tumor", drop = FALSE])
mean_normal <- rowMeans(log_expr[, condition == "Normal", drop = FALSE])

log2FC <- mean_tumor - mean_normal  # Tumor - Normal

# Adjusting p-values (FDR)
padj <- p.adjust(t_test_results, method = "BH")

# Combining into a results data frame
results_df <- data.frame(
  GeneID     = gene_info$GeneID,
  GeneSymbol = gene_info$GeneSymbol,
  log2FC     = log2FC,
  pvalue     = t_test_results,
  padj       = padj,
  stringsAsFactors = FALSE
)

# sort by adjusted p-value
results_df <- results_df %>% arrange(padj)

cat("\nTop 10 genes:\n")
print(head(results_df, 10))


# Save results to CSV

if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)

write.csv(results_df,
          file = "data/processed/breast_cancer_DE_results_simple_ttest.csv",
          row.names = FALSE)

cat("\nDifferential expression results saved to data/processed/breast_cancer_DE_results_simple_ttest.csv\n")


# Volcano plot

volcano_df <- results_df %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-10),  # avoid log(0)
    Significant = ifelse(padj < 0.05 & abs(log2FC) > 1, "Yes", "No")
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_padj, color = Significant)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  xlab("log2 Fold Change (Tumor - Normal)") +
  ylab("-log10(adjusted p-value)") +
  ggtitle("Volcano Plot: Tumor vs Normal")

print(p_volcano)

ggsave("figures/volcano_plot_tumor_vs_normal.png", p_volcano,
       width = 6, height = 4, dpi = 300)


# PCA plot of samples

# PCA on log-expression (transpose because prcomp expects samples in rows)
pca <- prcomp(t(log_expr), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  SampleID = colnames(log_expr),
  Condition = condition,
  Subtype   = meta_df$Subtype
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                            color = Condition,
                            shape = Subtype)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Samples (log2 expression)") +
  xlab("PC1") +
  ylab("PC2")

print(p_pca)

ggsave("figures/pca_tumor_normal.png", p_pca,
       width = 6, height = 4, dpi = 300)


#Heatmap of top 50 genes

top50_genes <- results_df$GeneID[1:50]

# Subsetting log-expression
top50_idx <- match(top50_genes, gene_info$GeneID)
heat_mat <- log_expr[top50_idx, ]

# Adding rownames (gene symbols)
rownames(heat_mat) <- results_df$GeneSymbol[1:50]

# Annotating samples
annotation_col <- data.frame(
  Condition = condition,
  Subtype   = meta_df$Subtype
)
rownames(annotation_col) <- meta_df$SampleID

pheatmap(heat_mat,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes")

# Save heatmap as file (pheatmap draws directly)
png("figures/heatmap_top50_genes.png", width = 800, height = 800)
pheatmap(heat_mat,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes")
dev.off()

cat("\nAll plots saved in 'figures/' and results in 'data/processed/'.\n")

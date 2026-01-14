# debug_subtype_structure.R
suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

message("Querying subtypes...")
subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")

message("Subtype Dataframe Dimensions:")
print(dim(subtypes))

message("Subtype Column Names:")
print(colnames(subtypes))

message("First 5 rows of subtypes:")
print(head(subtypes, 5))

# Check specifically for PAM50-like columns
pam50_cols <- grep("PAM50", colnames(subtypes), value = TRUE, ignore.case = TRUE)
message("Columns comparing 'PAM50':")
print(pam50_cols)

if (length(pam50_cols) > 0) {
  message("Table of first PAM50 col:")
  print(table(subtypes[[pam50_cols[1]]]))
}

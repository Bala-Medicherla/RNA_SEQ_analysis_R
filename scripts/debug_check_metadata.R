# debug_check_metadata.R

if (!file.exists("data/processed/sample_metadata.rds")) {
  stop("Metadata file not found. Did step 00 finish?")
}

meta <- readRDS("data/processed/sample_metadata.rds")

print("Dimensions of metadata:")
print(dim(meta))

print("Column names:")
print(names(meta))

if ("BRCA_Subtype_PAM50" %in% names(meta)) {
  print("Table of BRCA_Subtype_PAM50:")
  print(table(meta$BRCA_Subtype_PAM50, useNA="always"))
} else {
  print("BRCA_Subtype_PAM50 column is MISSING.")
}

if ("group" %in% names(meta)) {
  print("Table of group:")
  print(table(meta$group, useNA="always"))
}
############################################################
# 05_session_info.R
#
# Goal:
#   Capture session information for reproducibility.
############################################################

dir.create("results", showWarnings = FALSE, recursive = TRUE)

session_path <- "results/sessionInfo.txt"
writeLines(capture.output(sessionInfo()), session_path)

message(paste("Session info saved to", session_path))
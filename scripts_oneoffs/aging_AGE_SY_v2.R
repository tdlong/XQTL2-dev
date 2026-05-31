BASE <- "/Users/tdlong/Desktop/age_May_2026"

SCAN_FILES <- c(
  file.path(BASE, "AGE_SY10_F.scan.txt"),
  file.path(BASE, "AGE_SY20_F.scan.txt"),
  file.path(BASE, "AGE_SY10_M.scan.txt"),
  file.path(BASE, "AGE_SY20_M.scan.txt")
)
SCAN_LABELS  <- c("SY10 Female", "SY20 Female", "SY10 Male", "SY20 Male")
SCAN_COLOURS <- c("pink", "red", "lightblue", "blue")
THRESHOLD    <- 6
OUT_FILE     <- file.path(BASE, "AGE_SY_v2_4scan.png")
FORMAT       <- "powerpoint"
OUT_HEIGHT_IN <- 5.8
PEAKS        <- NULL
GENES        <- NULL

source("scripts_oneoffs/plot_pseudoscan.R")

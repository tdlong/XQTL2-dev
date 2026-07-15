# aging_AGE_SY_v3.R — 5-panel Wald figure, the four AGE_SY v3 scans overlaid.
#
# SY10/SY20 x Female/Male on one 5-panel (per-chromosome) figure, Mb on the X.
# Males blue (SY20) / light blue (SY10); females red (SY20) / pink (SY10).
#
# Run from the XQTL2-dev repo ROOT, on the cluster (where process/AGE_SY_v3 is):
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_SY/common/plotting/aging_AGE_SY_v3.R
#
# Output: figures/AGE_SY_v3_4scan.png

BASE <- "process/AGE_SY_v3"

SCAN_FILES <- c(
  file.path(BASE, "AGE_SY10_F", "AGE_SY10_F.scan.txt"),
  file.path(BASE, "AGE_SY20_F", "AGE_SY20_F.scan.txt"),
  file.path(BASE, "AGE_SY10_M", "AGE_SY10_M.scan.txt"),
  file.path(BASE, "AGE_SY20_M", "AGE_SY20_M.scan.txt")
)
SCAN_LABELS   <- c("SY10 Female", "SY20 Female", "SY10 Male", "SY20 Male")
SCAN_COLOURS  <- c("pink", "red", "lightblue", "blue")
THRESHOLD     <- 6
FORMAT        <- "powerpoint"
OUT_HEIGHT_IN <- 5.8
PEAKS         <- NULL
GENES         <- NULL

dir.create("figures", showWarnings = FALSE)
OUT_FILE <- "figures/AGE_SY_v3_4scan.png"

source("scripts_oneoffs/AGE_SY/common/plotting/plot_pseudoscan.R")

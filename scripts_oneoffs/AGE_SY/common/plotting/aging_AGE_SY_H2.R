# aging_AGE_SY_H2.R — 5-panel Cutler H² figure, the four AGE_SY scans overlaid.
#
# Same figure as aging_AGE_SY_v3.R, but the y-axis is the Cutl_H2 column of the
# scan files instead of Wald_log10p. Same four scans, same colours: males blue
# (SY20) / light blue (SY10); females red (SY20) / pink (SY10).
#
# The y-axis is left free — no 0 floor, no threshold line — so windows where
# Cutler H² comes out negative are shown as they are.
#
# Run from the XQTL2-dev repo ROOT, on the cluster (where process/AGE_SY is):
#   module load R/4.2.2
#   Rscript scripts_oneoffs/AGE_SY/common/plotting/aging_AGE_SY_H2.R
# or submit it as a job with round3_v3_R12/plot_4scan_H2.sh.
#
# Output: figures/AGE_SY_4scan_CutlH2.png

BASE <- "process/AGE_SY"   # to plot a different scan folder, edit this

SCAN_FILES <- c(
  file.path(BASE, "AGE_SY10_F", "AGE_SY10_F.scan.txt"),
  file.path(BASE, "AGE_SY20_F", "AGE_SY20_F.scan.txt"),
  file.path(BASE, "AGE_SY10_M", "AGE_SY10_M.scan.txt"),
  file.path(BASE, "AGE_SY20_M", "AGE_SY20_M.scan.txt")
)
SCAN_LABELS   <- c("SY10 Female", "SY20 Female", "SY10 Male", "SY20 Male")
SCAN_COLOURS  <- c("pink", "red", "lightblue", "blue")
THRESHOLD     <- NULL
FORMAT        <- "powerpoint"
OUT_HEIGHT_IN <- 5.8
PEAKS         <- NULL
GENES         <- NULL

Y_COL         <- "Cutl_H2"
Y_LAB         <- expression("Cutler" ~ H^2)
Y_FLOOR_ZERO  <- FALSE

dir.create("figures", showWarnings = FALSE)
OUT_FILE <- file.path("figures", paste0(basename(BASE), "_4scan_CutlH2.png"))

source("scripts_oneoffs/AGE_SY/common/plotting/plot_pseudoscan.R")

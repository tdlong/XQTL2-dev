# make_splithalf_designs.R — split each AGE_SY design into odd and even replicates.
#
# Takes the four live design files and writes eight, one per treatment x half,
# into helpfiles/AGE_SY/splithalf/. Rows are dropped, nothing else is changed:
# REP keeps its real value (R1 stays 1, R3 stays 3), which is only safe because
# XQTL2 #32 fixed Heritability() to match replicates on label rather than on
# position. Against the pre-#32 pipeline this would silently use half the
# replicates with mismatched selection proportions.
#
# Odd/even rather than 1-6/7-12 because R1-R12 arrived over months -- a block
# split would put batch straight into the split-half error term.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/splithalf/make_splithalf_designs.R

SRC_DIR <- "helpfiles/AGE_SY"
OUT_DIR <- "helpfiles/AGE_SY/splithalf"
SCANS   <- c("AGE_SY10_F", "AGE_SY20_F", "AGE_SY10_M", "AGE_SY20_M")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (scan in SCANS) {
  src <- file.path(SRC_DIR, paste0(scan, ".test.txt"))
  if (!file.exists(src)) stop("design not found: ", src)
  d <- read.table(src, header = TRUE)

  for (half in c("odd", "even")) {
    keep <- if (half == "odd") d$REP %% 2 == 1 else d$REP %% 2 == 0
    sub  <- d[keep, , drop = FALSE]
    rownames(sub) <- NULL          # row names must be unique; renumber like the source

    # Both arms must survive the split, or the scan has nothing to contrast.
    nC <- sum(sub$TRT == "C"); nZ <- sum(sub$TRT == "Z")
    if (nC == 0 || nZ == 0)
      stop(scan, " ", half, ": one arm is empty (C=", nC, ", Z=", nZ, ")")
    if (!setequal(sub$REP[sub$TRT == "C"], sub$REP[sub$TRT == "Z"]))
      stop(scan, " ", half, ": C and Z arms hold different replicates")

    out <- file.path(OUT_DIR, paste0(scan, ".", half, ".txt"))
    write.table(sub, out)
    cat(sprintf("%-28s C=%2d Z=%2d  REP=%s\n", basename(out), nC, nZ,
                paste(sort(unique(sub$REP)), collapse = ",")))
  }
}
cat("\nWrote", length(SCANS) * 2, "design files to", OUT_DIR, "\n")

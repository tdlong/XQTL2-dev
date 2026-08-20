# make_designs.R — AGE_SY without replicates 8 and 9.
#
# Replicates 8 and 9 came from the May 2023 cage; every other AGE_SY replicate
# came from the November 2023 cage (helpfiles/AGE_2024/population_assignment.txt).
# Dropping them leaves a single-cage experiment, so cage can no longer contribute
# to anything downstream.
#
# It also leaves 10 replicates -- 1,2,3,4,5,6,7,10,11,12 -- which split evenly
# into odd (1,3,5,7,11) and even (2,4,6,10,12). The split-half error term stays
# balanced at 5 and 5, which a contiguous cut would not manage.
#
# Writes 12 designs per run: for each of the four treatments, the 10-replicate
# design and its two halves.
#
# Run from the XQTL2-dev repo ROOT:
#   Rscript scripts_oneoffs/AGE_SY/nov_only/make_designs.R

SRC_DIR <- "helpfiles/AGE_SY"
OUT_DIR <- "helpfiles/AGE_SY/nov_only"
SCANS   <- c("AGE_SY10_F", "AGE_SY20_F", "AGE_SY10_M", "AGE_SY20_M")
DROP    <- c(8, 9)                       # the May 2023 cage

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# In the female design files the Sex column is entirely "F", which read.table
# converts to logical FALSE and would write back as FALSE. Harmless to the scan,
# which never reads Sex, but it corrupts the file for anything that does.
read_design <- function(f) {
  d <- read.table(f, header = TRUE)
  for (v in names(d)) if (is.logical(d[[v]])) d[[v]] <- ifelse(d[[v]], "T", "F")
  d
}

emit <- function(d, out, label) {
  rownames(d) <- NULL
  nC <- sum(d$TRT == "C"); nZ <- sum(d$TRT == "Z")
  if (nC == 0 || nZ == 0)
    stop(label, ": one arm is empty (C=", nC, ", Z=", nZ, ")")
  if (!setequal(d$REP[d$TRT == "C"], d$REP[d$TRT == "Z"]))
    stop(label, ": C and Z arms hold different replicates")
  write.table(d, out)
  cat(sprintf("  %-34s C=%2d Z=%2d  REP=%s\n", basename(out), nC, nZ,
              paste(sort(unique(d$REP)), collapse = ",")))
}

for (scan in SCANS) {
  src <- file.path(SRC_DIR, paste0(scan, ".test.txt"))
  if (!file.exists(src)) stop("design not found: ", src)
  d <- read_design(src)
  keep <- d[!d$REP %in% DROP, , drop = FALSE]

  emit(keep, file.path(OUT_DIR, paste0(scan, ".no89.txt")), scan)
  for (half in c("odd", "even")) {
    sel <- if (half == "odd") keep$REP %% 2 == 1 else keep$REP %% 2 == 0
    emit(keep[sel, , drop = FALSE],
         file.path(OUT_DIR, paste0(scan, ".no89.", half, ".txt")),
         paste(scan, half))
  }
}
cat("\nWrote", length(SCANS) * 3, "design files to", OUT_DIR, "\n")

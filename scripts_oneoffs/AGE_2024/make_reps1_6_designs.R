# make_reps1_6_designs.R — cut the four AGE_SY designs down to replicates 1-6.
#
# The Aug 2024 pilot has 6 cages, so the like-for-like AGE_SY comparison is its
# first 6 replicates rather than all 12. Subsetting the design file is all that
# is needed -- the haplotypes are unchanged, and since XQTL2 #32 the pipeline
# takes REP from the design rather than as a positional index, so dropping rows
# is safe.
#
# Writes helpfiles/AGE_SY/<scan>_R1to6.test.txt next to the originals.
# Called by run_all.sh; can be run on its own from the repo root.

scans <- c("AGE_SY10_F", "AGE_SY20_F", "AGE_SY10_M", "AGE_SY20_M")
MAXREP <- 6

for (s in scans) {
  src <- file.path("helpfiles/AGE_SY", paste0(s, ".test.txt"))
  dst <- file.path("helpfiles/AGE_SY", paste0(s, "_R1to6.test.txt"))
  if (!file.exists(src)) stop("missing ", src)

  d <- read.table(src, header = TRUE)

  # In the female files the Sex column is all "F", which read.table converts to
  # logical FALSE and would write back as FALSE -- silently corrupting the
  # design. Any logical column here can only have come from an all-"F" or
  # all-"T" text column, so map it straight back. (colClasses is not the fix:
  # it is positional, and the leading row-name column shifts it by one.)
  for (v in names(d)) if (is.logical(d[[v]]))
    d[[v]] <- ifelse(d[[v]], "T", "F")

  keep <- d[d$REP <= MAXREP, ]

  # both arms have to survive, or the scan has nothing to contrast
  n_sel <- sum(keep$TRT == "Z"); n_con <- sum(keep$TRT == "C")
  if (n_sel == 0 || n_con == 0)
    stop(s, ": kept ", n_sel, " selected and ", n_con, " control rows")

  write.table(keep, dst)
  cat(sprintf("   %-12s %2d rows -> %2d  (%d selected, %d control, REP %s)\n",
              s, nrow(d), nrow(keep), n_sel, n_con,
              paste(sort(unique(keep$REP)), collapse = ",")))
}

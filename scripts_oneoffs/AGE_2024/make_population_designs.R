# make_population_designs.R — split each female design by source population.
#
# Both experiments used two cages built from the same eight founders, set up
# May 2023 and November 2023 (helpfiles/AGE_2024/population_assignment.txt).
# This writes one design per food x population so the Wald profiles can be
# compared, which the existing scans cannot do -- AGE_SY replicates 1-6 are
# entirely Nov, so the pilot-vs-AGE_SY figure has population confounded with food.
#
# Females only. Six designs:
#   lab_May   (4 reps)    lab_Nov   (2 reps)
#   SY10_F_May (2)        SY10_F_Nov (10)
#   SY20_F_May (2)        SY20_F_Nov (10)
#
# Called by run_population_scans.sh; runs on its own from the repo root.

MAY_SY   <- c(8, 9)          # AGE_SY replicates from the May cage
MAY_PILOT <- 1:4             # pilot replicates from the May cage

# In the female design files the Sex column is entirely "F", which read.table
# turns into logical FALSE and would write back as FALSE. Map logicals to text.
read_design <- function(f) {
  d <- read.table(f, header = TRUE)
  for (v in names(d)) if (is.logical(d[[v]])) d[[v]] <- ifelse(d[[v]], "T", "F")
  d
}

emit <- function(d, keep, out, label) {
  k <- d[d$REP %in% keep, ]
  n_sel <- sum(k$TRT == "Z"); n_con <- sum(k$TRT == "C")
  if (n_sel < 2 || n_con < 2)
    stop(label, ": only ", n_sel, " selected and ", n_con, " control rows")
  write.table(k, out)
  cat(sprintf("   %-18s %2d reps  (REP %s)\n", label, n_sel,
              paste(sort(unique(k$REP)), collapse = ",")))
}

cat("pilot (lab food, female):\n")
p <- read_design("helpfiles/AGE_Aug13_24/Ageing_Aug13.txt")
emit(p, MAY_PILOT, "helpfiles/AGE_Aug13_24/Ageing_Aug13.May.txt", "lab_May")
emit(p, setdiff(unique(p$REP), MAY_PILOT),
     "helpfiles/AGE_Aug13_24/Ageing_Aug13.Nov.txt", "lab_Nov")

cat("AGE_SY females:\n")
for (s in c("AGE_SY10_F", "AGE_SY20_F")) {
  d <- read_design(file.path("helpfiles/AGE_SY", paste0(s, ".test.txt")))
  emit(d, MAY_SY, file.path("helpfiles/AGE_SY", paste0(s, "_May.test.txt")),
       paste0(s, "_May"))
  emit(d, setdiff(unique(d$REP), MAY_SY),
       file.path("helpfiles/AGE_SY", paste0(s, "_Nov.test.txt")),
       paste0(s, "_Nov"))
}

cat("\nCaution when reading the result:\n")
cat("  - the AGE_SY May scans rest on 2 replicates, so 1 df for the between-\n")
cat("    replicate variance; only a systematic shift will show.\n")
cat("  - within the pilot, May was the WEAKER selection (median 0.084 against\n")
cat("    0.049 for Nov), so its May-vs-Nov contrast is not population alone.\n")

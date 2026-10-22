# loo_genome.R — leave-one-cage-out prediction, genome-wide.
#
# RUN ON HPC3, from /dfs7/adl/tdlong/fly_pool/XQTL2-dev:
#   sbatch temp_aging/analysis/loo_genome.sh
#
# The two outputs are small and land in temp_aging/numbers/, which is tracked.
# loo_genome.sh commits and pushes them itself when the job finishes, so there is
# nothing to stage by hand -- just sbatch here, then git pull on the laptop.
#
# WHY THIS NEEDS THE CLUSTER: the meansBySample files are ~257 MB each and there
# are four. Everything here is a function of the founder frequencies they hold,
# so the subsetting happens there and two small files come back through git.
#
# NO NEW SCANS ARE NEEDED. meansBySample already carries every replicate, and the
# fold logic is arithmetic on those frequencies. The pipeline Wald is not used:
# positions are ranked by h2 computed on the NINE TRAINING CAGES, which is both
# the quantity being aggregated and free of the held-out cage.
#
# TILING. By default tiles are 1 cM. Set TILE_MB=1 to tile in physical distance
# instead, 1 Mb, writing to *_mb.txt.gz. Two reasons to want both:
#
#   Invariance. These estimates are meant to be properties of the data, not of
#   how markers were sampled. If the cM and Mb tilings disagree, something is
#   wrong with the aggregation rather than with the genome.
#
#   Identification. In cM every interior tile has the same neighbours at the
#   same correlations, so the LD score is constant -- 89% of tiles sit within
#   10% of the mean. That removes exactly the variation LD score regression uses
#   to identify h2 from a slope, leaving only a ratio. Tiled in Mb the LD score
#   varies with local recombination rate, since the correlation between two
#   tiles is set by the cM between them, and the regression becomes available.
#   Both tilings therefore carry cM per tile, not just the tile index.
#
# WHAT EACH FOLD DOES
#   train = 9 cages, test = the 10th.
#   weights   e_f = mean over train of (Z - C)/i      -- per founder, per step
#   positions ranked by train h2, greedily clumped at CLUMP_CM
#   prediction at a chosen position: Dp_hat = i_test * e_f   (i_test is KNOWN,
#     so nothing is fitted to the held-out cage)
#   h2 three ways, as in LOO_fig panel B:
#     nine    = 100k sum_f e_tr^2 / p      (true + inflation)
#     cross   = 100k sum_f e_tr e_te / p   (unbiased: independent noise)
#     single  = 100k sum_f e_te^2 / p      (one cage, maximally inflated)
#
# THE TOTAL. Genome-wide integral of the per-tile cross-product h2, divided by
# the footprint width. The cross-product version is used so the integral does not
# inherit the squaring bias that produces the null-window floor -- which is the
# thing that made the floor correction ambiguous when computed within-sample.

suppressMessages(library(tidyverse))

SCANDIR  <- "process/AGE_SY/Scans"
CMSRC    <- "process/AGE_SY/AGE_SY_4scan_no89.txt.gz"   # supplies cM per step
DESDIR   <- "helpfiles/AGE_SY/nov_only"
TILE_MB  <- Sys.getenv("TILE_MB", "0") == "1"   # physical tiles instead of genetic
SUF      <- if (TILE_MB) "_mb" else ""
tile_of  <- function(pos, cM) if (TILE_MB) floor(pos / 1e6) else floor(cM)
OUT_POS  <- sprintf("temp_aging/numbers/loo_positions%s.txt.gz", SUF)   # panel A + per-fold h2
OUT_TILE <- sprintf("temp_aging/numbers/loo_tiles%s.txt.gz", SUF)       # per-tile h2, for the integral
OUT_U    <- sprintf("temp_aging/numbers/loo_u%s.txt.gz", SUF)           # per-tile per-founder u, for the joint fit
OUT_R    <- sprintf("temp_aging/numbers/loo_pos_reps%s.txt.gz", SUF)    # EVERY replicate at EVERY chosen position
DROP_REP <- c(8, 9)
NMAX     <- 60        # positions kept per fold; curves run to 30, 60 leaves room
CLUMP_CM <- 5         # greedy clumping window; r^2 < 0.1 at 4 cM in this population
AF_FLOOR <- 0.01      # as in heritability_rep(), keeps 1/p from exploding

cells <- expand_grid(sugar = c("SY10","SY20"), sex = c("F","M")) %>%
  mutate(scan = paste0("AGE_", sugar, "_", sex),
         file = file.path(SCANDIR, scan, paste0(scan, ".meansBySample.txt")),
         des  = file.path(DESDIR, sprintf("AGE_%s_%s.no89.txt", sugar, sex)))
miss <- c(cells$file[!file.exists(cells$file)], cells$des[!file.exists(cells$des)])
if (length(miss)) stop("missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

cm <- read.table(CMSRC, header = TRUE, sep = "\t") %>% as_tibble() %>%
  distinct(chr, pos, cM)
cat("cM map:", nrow(cm), "steps   tiling:", if (TILE_MB) "1 Mb" else "1 cM", "\n")

# greedy clumping, exactly as a GWAS clumps: take the best, drop everything
# within CLUMP_CM on that arm, repeat
clump <- function(z, win, n) {
  out <- integer(0); z <- z[order(-z$h2), ]
  while (nrow(z) && length(out) < n) {
    out <- c(out, z$id[1])
    z <- z[!(z$chr == z$chr[1] & abs(z$cM - z$cM[1]) < win), ]
  }
  out
}

POS <- list(); TILE <- list(); UTILE <- list(); REPS <- list()
for (r in seq_len(nrow(cells))) {
  cl <- cells[r, ]
  cat("\n==", cl$scan, "==\n")
  prop <- read.table(cl$des, header = TRUE) %>%
    filter(TRT == "Z") %>%
    transmute(REP = as.integer(REP), P = Proportion, i = dnorm(qnorm(1 - P)) / P)

  d <- read.table(cl$file, header = TRUE) %>% as_tibble() %>%
    filter(!REP %in% DROP_REP, TRT %in% c("C", "Z")) %>%
    pivot_wider(names_from = TRT, values_from = freq) %>%
    inner_join(prop, by = "REP") %>%
    inner_join(cm, by = c("chr", "pos")) %>%
    mutate(D = Z - C, e = D / i,
           k = if (cl$sex == "M") ifelse(chr == "chrX", 1, 2) else 2)
  reps <- sort(unique(d$REP))
  cat("  steps:", n_distinct(paste(d$chr, d$pos)), "  reps:", length(reps), "\n")
  chosen_pos <- list()     # every position any fold selects, for the per-replicate dump

  for (ho in reps) {
    tr <- d %>% filter(REP != ho) %>%
      group_by(chr, pos, cM, founder, k) %>%
      summarise(e_tr = mean(e), p = mean(C), .groups = "drop")
    te <- d %>% filter(REP == ho) %>%
      select(chr, pos, founder, e_te = e, D_te = D, i_te = i)
    j <- inner_join(tr, te, by = c("chr", "pos", "founder")) %>%
      mutate(pd = pmax(p, AF_FLOOR))

    # per-step h2, three ways
    st <- j %>% group_by(chr, pos, cM) %>%
      summarise(nine   = sum(100 * k * e_tr^2 / pd),
                cross  = sum(100 * k * e_tr * e_te / pd),
                single = sum(100 * k * e_te^2 / pd), .groups = "drop")

    # per 1 cM tile, for the genome-wide integral
    TILE[[length(TILE)+1]] <- st %>%
      group_by(chr, tile = tile_of(pos, cM)) %>%
      summarise(across(c(nine, cross, single), mean), cM = mean(cM),
                pos = mean(pos), .groups = "drop") %>%
      mutate(sugar = cl$sugar, sex = cl$sex, holdout = ho)

    # Per-tile, per-founder standardised effects, for the joint (COJO / GREML)
    # estimate.  u = e/sqrt(p), so u_tr * u_te is exactly the integrand summed
    # into `cross` above -- kept unsummed because the joint estimate needs each
    # founder's PROFILE along the chromosome, not the per-tile total.
    UTILE[[length(UTILE)+1]] <- j %>%
      group_by(chr, tile = tile_of(pos, cM), founder, k) %>%
      summarise(u_tr = mean(e_tr / sqrt(pd)), u_te = mean(e_te / sqrt(pd)),
                cM = mean(cM), pos = mean(pos), .groups = "drop") %>%
      mutate(sugar = cl$sugar, sex = cl$sex, holdout = ho)

    # positions: ranked on the TRAINING cages only, then clumped
    sel <- st %>% mutate(id = row_number(), h2 = nine)
    keep <- clump(sel, CLUMP_CM, NMAX)
    chosen <- sel[keep, ] %>% mutate(rank = seq_along(keep))

    POS[[length(POS)+1]] <- j %>%
      inner_join(chosen %>% select(chr, pos, rank), by = c("chr","pos")) %>%
      transmute(sugar = cl$sugar, sex = cl$sex, holdout = ho, rank, chr, pos, cM,
                founder, p,
                pred = i_te * e_tr, obs = D_te,          # panel A
                h2_nine   = 100 * k * e_tr^2 / pd,       # panel B, per founder
                h2_cross  = 100 * k * e_tr * e_te / pd,
                h2_single = 100 * k * e_te^2 / pd)
    chosen_pos[[length(chosen_pos)+1]] <- chosen %>% select(chr, pos)
    cat("    holdout", ho, ": kept", length(keep), "positions\n")
  }

  # Per-replicate values at every position any fold chose. POS above stores only
  # the MEAN over the nine, which is enough for the prediction itself but not for
  # anything needing the spread among replicates -- a variance, a jackknife, the
  # s^2/R correction applied to the nine. Those then need another cluster run,
  # which is what this avoids. One row per replicate x founder x position.
  REPS[[length(REPS)+1]] <- d %>%
    inner_join(bind_rows(chosen_pos) %>% distinct(chr, pos), by = c("chr", "pos")) %>%
    transmute(sugar = cl$sugar, sex = cl$sex, REP, chr, pos, cM, founder, k,
              p = C, z = Z, dp = D, e, i)
  cat("  per-replicate dump:", n_distinct(paste(bind_rows(chosen_pos)$chr, bind_rows(chosen_pos)$pos)),
      "distinct positions over all folds\n")
  rm(d); gc()
}

dir.create(dirname(OUT_POS), showWarnings = FALSE, recursive = TRUE)
write.table(bind_rows(POS),  gzfile(OUT_POS),  row.names=FALSE, quote=FALSE, sep="\t")
write.table(bind_rows(TILE), gzfile(OUT_TILE), row.names=FALSE, quote=FALSE, sep="\t")
write.table(bind_rows(UTILE), gzfile(OUT_U),    row.names=FALSE, quote=FALSE, sep="\t")
write.table(bind_rows(REPS), gzfile(OUT_R),    row.names=FALSE, quote=FALSE, sep="\t")
cat("\nWrote:", OUT_POS, ",", OUT_TILE, ",", OUT_U, "and", OUT_R, "\n")
cat("loo_genome.sh commits these and tries to push.\n")

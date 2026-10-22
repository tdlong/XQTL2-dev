# loo_curves.R — the two curves of Figure Yb, as a table in numbers/.
#
#   bash temp_aging/analysis/run_numbers.sh loo_curves
#
# CEILING, how much heritability the top N regions hold. From loo_split.txt.gz:
# five replicates rank the positions, four others estimate the response, the
# held-out replicate scores it, five random splits per fold. Nothing that chose
# a position contributes to its effect, so the cross product e_est x e_te has
# expectation t^2 and carries no winner's curse. The nine-replicate cross
# product does: ranking and effect come from the same nine, so its expectation
# is (t + bias) x t, about 10% high here.
#
# PREDICTOR, what an N-marker predictor explains in a held-out replicate. From
# loo_positions.txt.gz, the real predictor: effects from all nine. A test-set
# quantity, so the bias in the effects degrades it rather than inflating it
# (simulation 2026-10-22: 1.6% high, a small-sample term).
#
# The two use different rankings -- five replicates for the ceiling, nine for
# the predictor -- so the ceiling is conservative: with nine-replicate ranking
# those regions would hold slightly more.

suppressMessages(library(tidyverse))

POS <- "temp_aging/numbers/loo_positions.txt.gz"
SPL <- "temp_aging/numbers/loo_split.txt.gz"
YS  <- "temp_aging/numbers/yardstick.tsv"
TSV <- "temp_aging/numbers/loo_curves.tsv"
AF  <- 0.01   # as in loo_genome.R

P <- read.table(POS, header = TRUE, sep = "\t") %>% as_tibble() %>% mutate(trt = paste0(sugar, "_", sex))
S <- read.table(SPL, header = TRUE, sep = "\t") %>% as_tibble() %>% mutate(trt = paste0(sugar, "_", sex))
Y <- read.table(YS,  header = TRUE, sep = "\t") %>% as_tibble()
NMAX <- min(P %>% count(trt, holdout) %>% pull(n) %>% max() / n_distinct(P$founder),
            S %>% group_by(trt, holdout, split) %>% summarise(n = n_distinct(rank), .groups = "drop") %>% pull(n) %>% min())
NMAX <- floor(NMAX)
cat("curves run to N =", NMAX, "(5 cM clumping leaves this many positions)\n")
jk <- function(x) { n <- length(x); sqrt((n - 1) / n * sum((x - mean(x))^2)) }

# predictor, per held-out replicate and combined
folds <- map_dfr(1:NMAX, function(N) P %>% filter(rank <= N) %>% group_by(trt, holdout) %>%
  summarise(cr = sum(h2_cross), nn = sum(h2_nine), .groups = "drop") %>% mutate(N = N))
pred <- folds %>% group_by(trt, N) %>%
  summarise(predictor = mean(cr)^2 / mean(nn),
            predictor_se = jk(map_dbl(seq_len(n()), ~ mean(cr[-.x])^2 / mean(nn[-.x]))), .groups = "drop")

# ceiling, averaged over folds x splits
ceil <- map_dfr(1:NMAX, function(N) S %>% filter(rank <= N) %>% group_by(trt, holdout, split) %>%
  summarise(cross = sum(100 * k * e_est * e_te / pmax(p, AF)), .groups = "drop") %>%
  group_by(trt) %>% summarise(ceiling = mean(cross),
    ceiling_se = jk(map_dbl(unique(holdout), ~ mean(cross[holdout != .x]))), .groups = "drop") %>% mutate(N = N))

out <- pred %>% inner_join(ceil, by = c("trt", "N")) %>%
  left_join(Y %>% transmute(trt, trait = h2), by = "trt")
write.table(out, TSV, row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nheritability (% of phenotypic variance)\n\n")
print(as.data.frame(out %>% filter(N %in% c(10, 20, 30, NMAX)) %>%
  transmute(trt, N, `regions hold` = round(ceiling, 1), `predictor explains` = round(predictor, 1),
            `trait h2` = round(trait, 1), `hold/trait` = round(ceiling / trait, 2),
            `explains/hold` = round(predictor / ceiling, 2))), row.names = FALSE)
cat("\nWrote:", TSV, "\n")

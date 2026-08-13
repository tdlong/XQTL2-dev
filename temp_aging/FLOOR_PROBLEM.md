# The h² floor is sex-dependent, and `b` does not report it

Found by plotting h² against Wald per character (`h2_vs_wald.R`, figure
`h2_vs_wald.png`). Quantified in `floor_sensitivity.R`.

## What the plot shows

At Wald-null windows (max −log10 *P* < 2 over the four characters), where true
h² is ~0, the raw uncorrected h² is not the same in the two sexes:

| character | median raw h² at null windows | sd | reported *b* |
|---|---|---|---|
| SY10 female | 0.497 | 0.121 | 0.64 |
| SY20 female | 0.487 | 0.139 | 0.63 |
| SY10 male | 0.577 | 0.237 | 0.63 |
| SY20 male | 0.554 | 0.222 | 0.65 |

The male floor is ~0.07 higher and roughly twice as variable. The reported bias
*b* is the same in all four to two decimal places, so a floor built as one
function *F*(*b*) applied to all characters **cannot** remove this — it subtracts
the same amount from a male window and a female window that report the same *b*.

## Consequence

Three defensible ways to handle the floor give three answers:

| floor treatment | null median, F / M | sex share, all | sex share, autosomes |
|---|---|---|---|
| one global *F*(*b*) (as published) | −0.07 / +0.02 | 10.2% | 7.3% |
| global *F*(*b*), then recentre each character on its null median | 0 / 0 | 8.4% | 6.5% |
| separate isotonic *F*(*b*) per character | −0.014 / −0.069 | 8.7% | 7.3% |

The per-character isotonic fit centres the *mean* at zero, but male h² is
right-skewed, so it leaves the median negative — which is why it does not agree
with the median recentring. None of the three is obviously correct.

Per-chromosome male-vs-female comparisons move much more than the headline:

| arm | global *F* | median recentred | per-character *F* |
|---|---|---|---|
| chrX | M > F at 100% | 100% | 98% |
| chr3R | M > F at 29% | 13% | 9% |
| chr2L | M > F at 79% | 71% | 61% |
| chr2R | M > F at 67% | 49% | 24% |
| chr3L | M > F at 69% | 53% | 47% |

## What survives

- **chrX strongly male-biased** — holds under all three, and is expected from
  hemizygosity anyway.
- **chr3R female-biased** — holds under all three and strengthens under both
  corrections, because the residual bias runs against it.
- **chr2L male-biased** — direction holds but the magnitude is set by the floor
  treatment. Weak.
- **chr2R and chr3L** — do not survive. Both reverse direction between
  treatments. The earlier draft called 2R male-biased and 3L male-biased; neither
  claim is supportable.
- **Genome-wide sex share is 8–10%, autosomes 6.5–7.3%**, not 10.2% / 7.3% as a
  point value.

## The published interval understates this

The bootstrap in `partition_by_wald_rank.R` resamples null windows and refits the
floor, so it captures sampling noise *in* the floor. It always refits the same
global *F*(*b*) form, so it does not capture the choice of whether the floor may
differ by character. The [5.9, 15.9] interval on the top-20% sex+diet share is
therefore too narrow by an unknown amount.

## Worth reporting upstream

XQTL2 issue #34 added `Cutl_H2_bias`. The level is roughly right — global *F*
maps *b* onto 0.51–0.64 and the observed null h² is 0.49–0.58 — but *b* does not
track the difference between male and female scans, which is real and about 0.07
in h² units. Male scans have lower coverage per pool; if *b* is computed from
haplotype-frequency sampling variance it may not be picking that up.

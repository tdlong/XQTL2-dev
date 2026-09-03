# Changes to make in RESULTS.md

Written 3 Sep 2026, against the RESULTS text current on that date. Driven by
the METHODS revision of the same day. Nothing here has been applied to
RESULTS.md — the list is for Tony to work through.

## A. Required by the Methods change

The bootstrap over tiles was removed from the analysis. Resampling tiles
measured how much genomic locations differ from one another and reported it as
uncertainty about the value at a location. The error term is now the sex:diet
interaction pooled with the split-half replicate difference, five degrees of
freedom, entirely within a tile.

**A1 — paragraph 3, delete the three intervals.** Currently:

> Over the tiles with LWP exceeding 7.5, 13.2% (95% CI 10.0-19.1) of the
> variation among the four groups is sex-specific and 1.8% (1.1-3.6)
> diet-specific, with no detectable interaction (-0.4%, -0.8 to 0.1).

The point estimates 13.2% and 1.8% are unchanged and stand. The three
parenthetical intervals are gone. Replacement:

> Over the tiles with LWP exceeding 7.5, 13.2% of the variation among the four
> groups is sex-specific and 1.8% diet-specific. Against that error the sex
> component is present at 44% of those tiles individually and the diet
> component at 13%.

**A2 — paragraph 3, say what the error term is.** The sentence before it says
"after subtracting the pure error componet the subsetting allowed us to
estimate". It is no longer pure error alone. Suggested insertion after it:

> The sex by diet interaction is nowhere detectable, exceeding its 5% critical
> value at 6% of tiles, so it is pooled with the replicate difference into an
> error term on five degrees of freedom.

**A3 — paragraph 3, "sets of five" -> "two interleaved sets of five."**
Methods now says interleaved throughout; the odd/even split was a coding
artifact of how replicates happen to be numbered.

## B. Numbers to check — I could not reproduce these

**B1 — paragraph 2, the three heritability thresholds.** Text says a LWP of 2,
7.5 and 15 corresponds to heritabilities of 0.14%, 0.37% and 0.65%. Neither
current source gives those:

| source | LWP 2 | LWP 7.5 | LWP 15 |
|---|---|---|---|
| `h2_threshold.R`, current | 0.06 (floor, Wald<2) | 0.38 | 0.67 |
| fitted line of Supp Fig 1, `0.00789(T-7)` | 0.09 | 0.33 | 0.62 |
| RESULTS as written | 0.14 | 0.37 | 0.65 |

Two things going on. The 7.5 and 15 figures look like an older run of
`h2_threshold.R` — I changed that script on 3 Sep so a tile is read at its
max-Wald step rather than taking `max(H2)` and `max(Wald)` independently, which
moved LWP 15 from 0.68 to 0.67. That is a 0.01 shift and does not explain 0.65.
And `h2_threshold.R` emits no value at LWP 2 at all; its first column is the
floor over tiles whose Wald never reaches 2, which is a different quantity from
the fitted curve read off at 2. Worth settling which of the two sources the
paragraph means, and quoting that one consistently for all three.

## C. Verified, no change needed

Checked against the current data on 3 Sep:

- 268 tiles, 32% above LWP 7.5, 12% above 15 — correct (85/268, 32/268).
- X is 24% of the map, 3% of its tiles above 7.5, none above 15, against 41%
  and 16% for the autosomes — all correct.
- 83 of 203 autosomal tiles above 7.5 — correct.
- Males higher at 52% of those tiles — correct.
- 18 peaks, 80% interval on width 0.22 to 1.18 cM (140–775 kb) — correct.
- r = 0.98 between LWP and heritability — correct (0.979).

## D. Typos

- paragraph 3: "partititoned" -> "partitioned"
- paragraph 3: "componet" -> "component"
- paragraph 1: "an contemporaneus" -> "contemporaneous"

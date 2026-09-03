# Proposed changes to RESULTS.md

Against the version current 3 Sep 2026. **Nothing here is applied.** Mark each
one — approve / deny / modify — and I will make the ones you approve.

---

## A. Forced by the Methods revision

The tile bootstrap was removed from the analysis: resampling tiles measures how
much genomic locations differ from one another, and we were reporting that as
uncertainty about the value at a location. The error term is now the sex:diet
interaction pooled with the split-half replicate difference — five degrees of
freedom, entirely within a tile.

### A1. Paragraph 3 — delete the three intervals

**Current:** "Over the tiles with LWP exceeding 7.5, 13.2% (95% CI 10.0-19.1) of
the variation among the four groups is sex-specific and 1.8% (1.1-3.6)
diet-specific, with no detectable interaction (-0.4%, -0.8 to 0.1)."

**Proposed:** "Over the tiles with LWP exceeding 7.5, 13.2% of the variation
among the four groups is sex-specific and 1.8% diet-specific. Against that error
the sex component is present at 44% of those tiles individually and the diet
component at 13%."

The 13.2% and 1.8% are unchanged. Only the parentheticals go.

`[ ] approve   [ ] deny   [ ] modify:`

### A2. Paragraph 3 — say what the error term now is

**Insert** after "...after subtracting the pure error componet the subsetting
allowed us to estimate.":

"The sex by diet interaction is nowhere detectable, exceeding its 5% critical
value at 6% of tiles, so it is pooled with the replicate difference into an
error term on five degrees of freedom."

Also "pure error" in the preceding sentence is no longer accurate on its own.

`[ ] approve   [ ] deny   [ ] modify:`

### A3. Paragraph 3 — "sets of five" to "two interleaved sets of five"

Methods now says interleaved throughout. Odd/even was a coding artifact of how
the replicates happen to be numbered.

`[ ] approve   [ ] deny   [ ] modify:`

---

## B. Statements I could not support

### B1. Paragraph 2 — the three heritability thresholds

**Current:** "This function predicts heritabilities of 0.14%, 0.37% and 0.65% at
LWP thresholds of 2, 7.5 and 15 repectively."

Neither source produces those numbers:

| source | LWP 2 | LWP 7.5 | LWP 15 |
|---|---|---|---|
| `h2_threshold.R`, current | 0.06 (floor, Wald<2) | 0.38 | 0.67 |
| fitted line of Supp Fig 1, `0.00789(T-7)` | 0.09 | 0.33 | 0.62 |
| **RESULTS as written** | **0.14** | **0.37** | **0.65** |

Two separate problems. The 7.5 and 15 figures look like an older run of
`h2_threshold.R`; I changed that script on 3 Sep so a tile is read at its
max-Wald step rather than taking `max(H2)` and `max(Wald)` independently, which
moved LWP 15 from 0.68 to 0.67 — a 0.01 shift that does not explain 0.65. And
`h2_threshold.R` emits nothing at LWP 2: its first column is a floor over tiles
whose Wald never reaches 2, a different quantity from the fitted curve read off
at 2.

**Needs a decision on which source the sentence means**, then all three numbers
quoted from that one source. I can produce either set.

`[ ] use h2_threshold.R   [ ] use the fitted line   [ ] leave, I will check`

### B2. Paragraph 1 — the two width ranges are not the same interval

**Current:** "an 80% confidence interval on width of 0.22 to 1.18 cM (140 kb to
775 kb)."

Both pairs are correct as 10th/90th percentiles, but they are percentiles of two
different distributions and describe different peaks. The 0.22 cM peak is chr3L
9.30 Mb, which is 55 kb. The 140 kb peak is chr3R 14.64 Mb, which is 0.26 cM.
Written with the kb in parentheses it reads as one interval converted to other
units, which it is not.

**Proposed:** "...an 80% range on width of 0.22 to 1.18 cM, or equivalently 140
to 775 kb." — or drop the kb figures.

`[ ] approve   [ ] deny   [ ] modify:`

### B3. Paragraph 4 — "most tightly localised" holds in kb but not cM

**Current:** "It is the most tightly localised of the peaks in Table S1, with a
two LWP support interval of 55 kb or 0.22 cM."

55 kb is the narrowest in Table S1. But in cM it is not: chr3L 21.62 Mb has a
support interval of 0.16 cM against this peak's 0.22 cM. Since the paper argues
resolution is genetic, a reader may check the cM column.

**Proposed:** "It has the narrowest support interval of any peak in Table S1 in
physical terms, 55 kb, or 0.22 cM."

`[ ] approve   [ ] deny   [ ] modify:`

### B4. Paragraph 2 — the 2% claim is asserted, not calculated

**Current:** "Regions contributing less than 2% to heritable variation would
have gone undetected in prior work employing modestly sized mapping panels
consisting of several hundred RILs (cites)..."

Nothing in the analysis establishes the 2% figure or the RIL comparison. The
Methods now carry the machinery to compute it — the non-centrality relation
gives the heritability detectable at a given power for a given design — so this
could be made a real number for a few hundred RILs rather than an assertion.

`[ ] compute it   [ ] soften the wording   [ ] leave as is`

### B5. Paragraph 2 — where 0.4% comes from

**Current:** "Despite 0.4% being a subtle heritability estimate at the scale of
the experiment we carried out..."

0.4% does not match either figure quoted in the preceding sentence (0.37% at LWP
7.5). Probably meant as a round number for that, but it reads as a third value.
Resolves itself once B1 is settled.

`[ ] approve tying it to the B1 number   [ ] deny   [ ] modify:`

---

## C. Typos

- P1: "an contemporaneus" -> "contemporaneous"
- P2: "repectively" -> "respectively"
- P3: "partititoned" -> "partitioned"
- P3: "componet" -> "component"
- P4: "a much larger allele frequency change that is observed" -> "than is observed"
- P4: "A large-effect loci sits" -> "locus" (or "Large-effect loci sit")

`[ ] approve all   [ ] deny   [ ] modify:`

---

## D. Checked against the data, correct as written

No change needed. Verified 3 Sep 2026 against the current scan and split-half
files:

- 268 tiles; 32% above LWP 7.5 (85/268), 12% above 15 (32/268).
- Selected proportions 3–11% (measured 3.1% to 10.8%), mean ~5%.
- X is 24% of the map; 3% of its tiles above 7.5 (2/65); none above 15;
  autosomes 41% and 16%.
- 83 of 203 autosomal tiles above 7.5.
- Males higher at 52% of those tiles.
- 18 peaks in Table S1 (8 above LWP 15, 10 between 7.5 and 15).
- r = 0.98 between LWP and heritability (0.979).
- chr3L 9.30 Mb: LWP 163.8 in SY20 males, minimum 36.1 across the four groups;
  h2 4.69% SY20 male and 1.53% SY10 female; next largest peak 1.67%, so
  "two to three times" holds at 2.8x.

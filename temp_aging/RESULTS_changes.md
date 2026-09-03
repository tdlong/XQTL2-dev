# Proposed changes to RESULTS.md — round 2

Round 1 items A (all three) and C (typos) are **applied**. What follows is the
B set, recalculated as you asked. **None of B is applied.** Mark each and I
will make the ones you approve.

---

## Applied in this round, for the record

- **A1** the three bootstrap intervals deleted; 13.2% and 1.8% kept; added
  "Against that error the sex component is present at 44% of those tiles
  individually and the diet component at 13%."
- **A2** inserted "The sex by diet interaction is nowhere detectable, exceeding
  its 5% critical value at 6% of tiles, so it is pooled with the replicate
  difference into an error term on five degrees of freedom." and changed
  "pure error componet" to "error component".
- **A3** "sets of five" -> "two interleaved sets of five".
- **C** contemporaneus -> contemporaneous; repectively -> respectively;
  partititoned -> partitioned; componet -> component; "change that is observed
  at any other peaks" -> "change than is observed at any other peak";
  "A large-effect loci sits" -> "A large-effect locus sits".

---

## B1. The three heritability thresholds — recalculated

**Current:** "This function predicts heritabilities of 0.14%, 0.37% and 0.65% at
LWP thresholds of 2, 7.5 and 15 respectively."

You are right that when two sources disagree we recalculate. The sentence says
*this function predicts*, so the source has to be the function in the Methods —
`h2 = c(T - 7)`, not the empirical tile averages of `h2_threshold.R`, which
answer a different question (what heritability does a tile at that Wald carry on
average).

Refitted just now against the current scan, 89,680 step x treatment values:

| | c | LWP 2 | LWP 7.5 | LWP 15 |
|---|---|---|---|---|
| **refit, current data** | **0.00774** | **0.09%** | **0.32%** | **0.61%** |
| previous fit | 0.00789 | 0.09 | 0.33 | 0.62 |
| RESULTS as written | — | 0.14 | 0.37 | 0.65 |

**Proposed:** "...predicts heritabilities of 0.09%, 0.32% and 0.61% at LWP
thresholds of 2, 7.5 and 15 respectively."

`make_FigS1.R` refits this constant each run rather than hardcoding it, so the
figure already uses 0.00774 and needs no change.

`[ ] approve   [ ] deny   [ ] modify:`

## B2. The peak width range — recalculated

**Current:** "an 80% confidence interval on width of 0.22 to 1.18 cM (140 kb to
775 kb)."

Recalculated with `quantile()` on the 18 Table S1 peaks, which reproduces what
`peak_table.R` prints:

| | 10th | median | 90th |
|---|---|---|---|
| cM | **0.25** | 0.71 | **1.51** |
| kb | **168** | 252 | **828** |

So both pairs in the text are stale. My round-1 point about the two ranges was
a separate and smaller thing: they are quantiles of two different distributions,
so the 0.25 cM peak and the 168 kb peak are not the same peak. That is fine as
long as the sentence does not read as one interval converted into other units.

**Proposed:** "...results in an 80% range on width of 0.25 to 1.51 cM, and of
168 to 828 kb."

One loose end: `peak_table.R` labels this output "pooled over all 23 peaks"
while Table S1 and the RESULTS text both say 18. The numbers match the 18, so
the label looks stale, but worth confirming which is right.

`[ ] approve   [ ] deny   [ ] modify:`

## B3. "Most tightly localised" — the context you asked for

**Current:** "It is the most tightly localised of the peaks in Table S1, with a
two LWP support interval of 55 kb or 0.22 cM."

Table S1 gives every peak's support interval twice, once in kb and once in cM.
The chr3L 9.30 Mb peak is 55 kb, and 55 kb is the smallest number in the kb
column — so in physical distance it is indeed the most tightly localised.

But it is 0.22 cM, and the chr3L 21.62 Mb peak is 0.16 cM. In genetic distance
that peak is tighter, so "most tightly localised" is not true if the reader
reads the cM column. It matters because the paper argues elsewhere that mapping
resolution is genetic, not physical, which invites exactly that reading.

**Proposed:** "It has the narrowest support interval of any peak in Table S1 in
physical terms, 55 kb, or 0.22 cM."

`[ ] approve   [ ] deny   [ ] modify:`

## B4. The 2% claim — the context you asked for, and a number

**Current:** "Regions contributing less than 2% to heritable variation would
have gone undetected in prior work employing modestly sized mapping panels
consisting of several hundred RILs (cites)..."

The claim is that a study with a few hundred RILs could not have found what we
find. Nothing in our analysis produces the 2% figure — it is asserted. Since the
Methods now carry the power relation, it can be computed instead.

For a one-degree-of-freedom test at a LOD-3 threshold with 80% power, the
non-centrality needed is 22.4, and a QTL explaining a fraction h2 of line
variance in n RILs gives a non-centrality of n*h2. So the smallest detectable
effect is 22.4/n:

| RILs | smallest detectable |
|---|---|
| 200 | 11.2% |
| 300 | 7.5% |
| 500 | 4.5% |
| 1000 | 2.2% |

So 2% is much too generous: at several hundred RILs the real floor is roughly
5–11%, and 2% is only reached at about a thousand lines. The claim is right in
direction and understated in magnitude.

**Proposed:** "Regions contributing less than roughly 5% of variation would have
gone undetected in prior work employing modestly sized mapping panels consisting
of several hundred RILs (cites)..."

`[ ] approve 5%   [ ] use a different figure:   [ ] leave as 2%`

## B5. The 0.4%

You are right, it is a rounding of the LWP 7.5 value. If B1 is approved that
value becomes 0.32%, so this should follow it rather than stay at 0.4%.

**Proposed:** "Despite 0.32% being a subtle heritability estimate at the scale
of the experiment we carried out..."

`[ ] approve   [ ] deny   [ ] modify:`

---

## Checked and correct, no change

Verified 3 Sep 2026 against the current scan and split-half files: 268 tiles;
32% above LWP 7.5 and 12% above 15; selected proportions 3–11% (3.1 to 10.8),
mean ~5%; X is 24% of the map with 3% of its tiles above 7.5 and none above 15,
against 41% and 16% for the autosomes; 83 of 203 autosomal tiles above 7.5;
males higher at 52% of those; 18 peaks (8 above 15, 10 between 7.5 and 15);
r = 0.98 between LWP and heritability (0.979); chr3L 9.30 Mb reaching LWP 163.8
in SY20 males with a minimum of 36.1 across the four groups, h2 of 4.69% in SY20
males and 1.53% in SY10 females, and "two to three times the largest of the
other regions" holding at 2.8x against the next largest peak at 1.67%.

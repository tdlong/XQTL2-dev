# RESULTS.md changes — record

All items from rounds 1–3 are **applied**. One new item is open at the bottom.

## Applied

### Round 1 — forced by the Methods revision

- **A1** The three bootstrap intervals deleted from paragraph 3. The tile
  bootstrap was removed from the analysis: resampling tiles measured how much
  genomic locations differ from one another and reported it as uncertainty about
  the value at a location. The 13.2% and 1.8% point estimates stand; added
  "the sex component is present at 44% of those tiles individually and the diet
  component at 13%".
- **A2** Added the pooled error term — the sex by diet interaction pooled with
  the replicate difference, five degrees of freedom.
- **A3** "sets of five" → "two interleaved sets of five".

### Round 2 — statements that could not be supported

- **B1** Heritability thresholds recalculated from the refitted relation
  (`c = 0.00774` on the current scan, 89,680 step × treatment values):
  **0.09%, 0.32%, 0.61%** at LWP 2, 7.5 and 15, replacing 0.14/0.37/0.65.
  `make_FigS1.R` refits rather than hardcodes, so the figure was already current.
- **B2** Peak widths recomputed with `quantile()` on the 18 Table S1 peaks:
  **0.25 to 1.51 cM** and **168 to 828 kb**, replacing 0.22–1.18 and 140–775.
  `peak_table.R` also had "23 peaks" and an ordinal for 23 hardcoded while
  pooling 18; both now derive from `nrow(tab)`.
- **B3** "most tightly localised" reworded to Tony's version — narrowest in
  physical terms, 55 kb, cM dropped. The claim held in kb but failed in cM,
  where chr3L 21.62 Mb is narrower at 0.16 cM.
- **B5** 0.4% → **0.32%**, following B1.

### Round 3 — the 2% claim

- **B4** The placeholder parenthetical ("a like-for-like power comparison is
  still to be done") replaced, since the comparison is now done:
  `design_power.R`, `DESIGN_POWER.md`, `make_FigS2.R` → `FigureS2_plot.png`.
  New text concedes the per-fly inefficiency, gives the mechanism, and points at
  the supplement. Tony has this to edit.

### Round 1 — typos

contemporaneus → contemporaneous; repectively → respectively; partititoned →
partitioned; componet → component; "change that is observed at any other peaks"
→ "than is observed at any other peak"; "A large-effect loci sits" → "locus".

## Open

### D1. Paragraph 3 repeats the loose interaction wording

Paragraph 3 currently says:

> The sex by diet interaction is nowhere detectable, exceeding its 5% critical
> value at 6% of tiles, so it is pooled with the replicate difference…

This is the phrasing Tony objected to in the Methods — 6% is not 5%, and the
sentence gives no reason to treat the difference as noise. The Methods was fixed
to lead with the distribution and give counts rather than a percentage. RESULTS
should match:

> The sex by diet interaction is nowhere detectable, five of the eighty-five
> tiles exceeding the 5% critical value where four are expected, so it is pooled
> with the replicate difference…

For the record: median `MS_sex:diet / MS_rep` over significant tiles is 0.46
against a null median of 0.55, and P(X ≥ 5 | n = 85, p = 0.05) = 0.42.

`[ ] approve   [ ] deny   [ ] modify:`

## Checked and correct, no change

Verified against the current scan and split-half files: 268 tiles; 32% above
LWP 7.5 and 12% above 15; selected proportions 3–11% (3.1 to 10.8), mean ~5%;
X is 24% of the map with 3% of its tiles above 7.5 and none above 15, against
41% and 16% for the autosomes; 83 of 203 autosomal tiles above 7.5; males higher
at 52% of those; 18 peaks (8 above 15, 10 between 7.5 and 15); r = 0.98 between
LWP and heritability (0.979); chr3L 9.30 Mb reaching LWP 163.8 in SY20 males
with a minimum of 36.1 across the four groups, h2 of 4.69% in SY20 males and
1.53% in SY10 females, and "two to three times the largest of the other regions"
holding at 2.8x against the next largest peak at 1.67%.

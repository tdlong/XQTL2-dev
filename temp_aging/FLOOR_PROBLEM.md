> **SUPERSEDED — do not quote.** Written 13 August 2026 against 12 replicates, 75 kb windows and a Wald threshold of 5. The analysis now uses 10 replicates (8 and 9 dropped), 268 one cM tiles, and thresholds of 7.5 and 15. See MANIFEST.md and results_para1-4.md.

# Note on the h² floor: where it matters and where it does not

*h*² squares founder effects estimated with error, so a window at which nothing
happened still returns a positive value — the floor, roughly 0.5 on the Fig. 1b
scale. `h2_threshold.R` and `significant_regions.R` between them settle how much
this matters, and the answer depends entirely on where you look.

## Where it does not matter

Restricted to the significant genome (Wald > 5 in at least one group), the
male-versus-female comparison is **unchanged** by floor removal:

| arm | M > F at, floor removed | floor left in |
|---|---|---|
| chrX | 100% | 100% |
| chr2L | 82% | 82% |
| chr3L | 81% | 81% |
| chr2R | 59% | 59% |
| chr3R | 19% | 19% |

The floor is near-constant across groups at a given window, so it cancels from
any within-window comparison between them. Every sex claim in the draft rests on
these regions and none of them depends on how the floor was estimated.

## Where it does matter

The partition into shared / sex / diet components does move — sex is 10.9% with
the floor removed and 3.8% with it left in — because the floor adds equally to
all four groups and so inflates the shared component. This is why Fig. 1c is
floor-corrected and Fig. 1b is not.

It matters more the less signal there is. Over the whole genome, including the
45% that never clears the threshold, the floor dominates and small changes in how
it is fitted swing the answer around; three defensible treatments give sex shares
of 10.2%, 8.4% and 8.7%. Restricted to significant tiles the estimate is stable
because the quantity being corrected is large relative to the correction.

## Worth reporting upstream

XQTL2 issue #34 added `Cutl_H2_bias`. Its level is about right — the global fit
maps *b* onto 0.51–0.64 and observed null h² is 0.49–0.58 — but *b* does not
distinguish the male from the female scans, reading 0.63–0.65 in all four while
the actual null distributions differ (95th percentile 1.06–1.10 in males against
0.70–0.74 in females). Male scans have lower coverage per pool; if *b* comes from
haplotype-frequency sampling variance it may not be capturing that. This does not
affect any conclusion here, for the reason above, but a user working closer to
the floor would be misled.

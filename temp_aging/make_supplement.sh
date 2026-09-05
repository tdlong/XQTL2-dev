#!/bin/bash
# make_supplement.sh -- build the supplementary .docx for the Longevity manuscript.
#
#   bash temp_aging/make_supplement.sh
#
# Writes temp_aging/Longevity_supplement.docx, for uploading to Drive and opening
# with Google Docs. Regenerate rather than editing the .docx or any intermediate.
#
# Column widths: pandoc derives relative widths for a docx table from the widths
# of the dashes in the markdown separator row, so the python below sizes those to
# the widest cell in each column. Section subtitles are paragraphs BETWEEN tables,
# never rows inside one -- a long string in a cell makes Word autofit blow the
# column out.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=temp_aging/Longevity_supplement.docx
MD=$(mktemp -t supp).md

python3 - "$MD" <<'PY'
import csv, sys, subprocess
md = open(sys.argv[1], 'w')
sha = subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip()
import datetime; today = datetime.date.today().strftime('%-d %b %Y')

W = md.write
W("# Longevity — Supplementary material\n\n")
W(f"From the XQTL2-dev analysis repo, commit {sha}, {today}. Figures are PNGs in the\n"
  "repo and must be inserted by hand; everything else here is final text.\n\n")

# ---------------- Table S1 ----------------
W("## Supplementary Table 1\n\n")
W("Peaks with a local maximum LWP exceeding 15 or 7.5, each at least 5 cM from any\n"
  "stronger peak on the same arm. The interval is the two-LWP support interval, walked\n"
  "out independently in each direction and so **not centred on the peak** — at 3L 21.62\n"
  "it runs 0.25 Mb to the left and 1.24 Mb to the right. Width is given as **kb / cM**\n"
  "and each treatment column as **-log10 P / h² (%)**. The treatment with the strongest\n"
  "signal is the largest of the four, so it has no column of its own.\n\n")

rows = list(csv.reader(open('temp_aging/TableS1.tsv'), delimiter='\t'))
HDR = ["Peak","cM","Interval (Mb)","Width","SY10 F","SY20 F","SY10 M","SY20 M"]
# every column left-aligned: a right-flush number beside a left-flush one
# collides at the column boundary however wide the columns are
RIGHT = set()
sections, cur = [], None
for r in rows[1:]:
    if not r or not r[0].strip(): continue
    if not r[0].startswith('chr'):
        cur = {"title": r[0], "rows": []}; sections.append(cur); continue
    cur["rows"].append([f"{r[0].replace('chr','')} {r[1]}", r[2], r[3],
                        f"{r[4]} / {r[5]}"] + [f"{r[7+i]} / {r[11+i]}" for i in range(4)])

allr = [x for s in sections for x in s["rows"]]
# +2 characters of slack per column: rendered width is not proportional to
# character count (digits are narrower than letters) and columns that abut with
# no gap read as one number running into the next.
w = [max([len(HDR[i])] + [len(r[i]) for r in allr]) + 2 for i in range(len(HDR))]
sep = "|" + "|".join(("-"*n+":") if i in RIGHT else (":"+"-"*n) for i,n in enumerate(w)) + "|"
row = lambda c: "| " + " | ".join(x.ljust(w[i]) for i,x in enumerate(c)) + " |"
for s in sections:
    W(f"**{s['title']}**\n\n{row(HDR)}\n{sep}\n")
    for r in s["rows"]: W(row(r) + "\n")
    W("\n")

# ---------------- figures ----------------
W("## Supplementary Figure 1\n\n*File: FigureS1_plot.png, from make_FigS1.R*\n\n")
W("Heritability against the Wald statistic, all four treatments pooled, with the\n"
  "relation h² = c(T − 7) drawn as a line. At non-significant steps heritability is\n"
  "over-estimated; over moderately significant LWPs of 5 to 15 it is over-estimated by\n"
  "an average of 0.05 percentage points, a 16–17% relative inflation; at highly\n"
  "significant steps the estimate is accurate.\n\n")
W("## Supplementary Figure 2\n\n*File: FigureS2_plot.png, from make_FigS2.R*\n\n")
W("Power against experiment size for four designs, at a locus explaining 2% of\n"
  "phenotypic variance (left) and 0.5%, the size resolved here (right). The x axis is\n"
  "individuals phenotyped, or lines for the inbred-line and RIL panels. Solid segments\n"
  "span the sizes these experiments are actually run at; dotted segments are\n"
  "hypothetical extensions that nobody has run. Points mark each design as currently\n"
  "practised, and the dashed horizontal line is 80% power.\n\n")

# ---------------- Supplementary Result X ----------------
txt = open('temp_aging/DESIGN_POWER.md').read()
txt = txt.replace("# What conventional designs could have detected",
                  "## Supplementary Result X — what conventional designs could have detected", 1)
W(txt)
md.close()
PY

pandoc "$MD" --reference-doc=temp_aging/reference.docx -o "$OUT"
python3 temp_aging/style_tables.py "$OUT"
rm -f "$MD"
echo "Wrote $OUT"

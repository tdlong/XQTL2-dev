#!/usr/bin/env python3
"""
make_AGE_SY_design_files.py — Generate AGE_SY pipeline design files from summary_info_v1.xlsx

Reads helpfiles/AGE_SY/summary_info_v1.xlsx (two sheets):

  Control_Flies:
    Replicate           — Rep1, Rep2, ... (skip RepA* rows — different base population)
    Sex                 — F or M
    Num_Random_Controls — fly count (used as Num in design file)

  Aged_Flies:
    Replicate           — Rep1, Rep2, ... (skip RepA* rows)
    Sex                 — F or M
    Diet                — SY10 or SY20
    Num_Aged            — fly count (used as Num in design file)
    Percent_Aged        — divide by 100 → Proportion in design file

BAM name convention: {longTRT}_R{n}_{sex}
  e.g. AgeSY10_R8_F, AgeSY20_R9_M, Con_R10_F

Output files (in helpfiles/AGE_SY/):
  AGE_SY10_F.test.txt  — AgeSY10 vs Con, females, R1-R12
  AGE_SY10_M.test.txt  — AgeSY10 vs Con, males,   R1-R12
  AGE_SY20_F.test.txt  — AgeSY20 vs Con, females, R1-R12
  AGE_SY20_M.test.txt  — AgeSY20 vs Con, males,   R1-R12

IMPORTANT: filesize (column 2) is written as 0 — a placeholder.
           After the BAMs are merged, update with actual sizes:
             stat -c "%s" data/bam/AGE_SY/<sample>.bam   (Linux)
           Then replace the 0s in each design file (or rerun this script
           with the --filesizes option once that is implemented).

Run from: XQTL2-dev root
  python3 scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py

NOTE: R12 (round 3, July_26 run) requires its fly counts to be present in
      summary_info_v1.xlsx first; rows missing from the xlsx are warned and
      skipped, so add R12 to the sheets before relying on the R12 output.
"""

import sys
import openpyxl
from pathlib import Path

XLSX   = Path("helpfiles/AGE_SY/summary_info_v1.xlsx")
OUTDIR = Path("helpfiles/AGE_SY")
REPS   = list(range(1, 13))   # R1–R12

# ── Read Control_Flies ────────────────────────────────────────────────────────
# control[(rep_num, sex)] = Num_Random_Controls
control = {}
ws = openpyxl.load_workbook(XLSX, data_only=True)["Control_Flies"]
first = True
for row in ws.iter_rows(values_only=True):
    if first:
        first = False
        continue
    rep_str, _, sex, _, num = row[:5]
    if not rep_str or not str(rep_str).startswith("Rep"):
        continue
    rep_raw = str(rep_str)[3:]          # "1", "2", ..., "11"
    if not rep_raw.isdigit():           # skip RepA1, RepA2, etc.
        continue
    control[(int(rep_raw), sex)] = num

# ── Read Aged_Flies ───────────────────────────────────────────────────────────
# aged[(rep_num, sex, diet)] = (Num_Aged, Proportion)
aged = {}
ws = openpyxl.load_workbook(XLSX, data_only=True)["Aged_Flies"]
first = True
for row in ws.iter_rows(values_only=True):
    if first:
        first = False
        continue
    rep_str, _, sex, diet, _, _, num_aged, pct_aged = row[:8]
    if not rep_str or not str(rep_str).startswith("Rep"):
        continue
    rep_raw = str(rep_str)[3:]
    if not rep_raw.isdigit():
        continue
    rep_num = int(rep_raw)
    if num_aged in (None, "-") or pct_aged in (None, "-"):
        aged[(rep_num, sex, diet)] = (None, None)
    else:
        aged[(rep_num, sex, diet)] = (int(num_aged), round(float(pct_aged) / 100, 4))

# ── Build rows for one design file ───────────────────────────────────────────
def make_rows(diet, sex):
    long_trt = f"AgeSY{diet[2:]}"      # "SY10" → "AgeSY10"
    rows = []

    # Treatment rows first, then control rows
    for long_trt_local, trt_code, reps_iter in [
        (long_trt, "Z", REPS),
        ("Con",    "C", REPS),
    ]:
        for rep in reps_iter:
            if trt_code == "Z":
                key = (rep, sex, diet)
                if key not in aged:
                    print(f"WARNING: no aged data for {key}", file=sys.stderr)
                    continue
                num, prop = aged[key]
                prop_str = "NA" if prop is None else f"{prop:.4f}"
                num_str  = "NA" if num  is None else str(num)
            else:
                key = (rep, sex)
                if key not in control:
                    print(f"WARNING: no control data for {key}", file=sys.stderr)
                    continue
                num_str  = str(control[key])
                prop_str = "NA"

            bam = f"{long_trt_local}_R{rep}_{sex}"
            rows.append({
                "bam":     bam,
                "longTRT": long_trt_local,
                "Rep":     f"R{rep}",
                "Sex":     sex,
                "TRT":     trt_code,
                "REP":     rep,
                "Num":     num_str,
                "Prop":    prop_str,
            })
    return rows

# ── Write one design file ─────────────────────────────────────────────────────
HEADER = '"filesize" "bam" "longTRT" "Replicate" "Sex" "TRT" "REP" "REPrep" "Num" "Proportion"'

def write_design(path, rows):
    lines = [HEADER]
    for i, r in enumerate(rows, 1):
        lines.append(
            f'"{i}" 0 "{r["bam"]}" "{r["longTRT"]}" "{r["Rep"]}" '
            f'"{r["Sex"]}" "{r["TRT"]}" {r["REP"]} 1 {r["Num"]} {r["Prop"]}'
        )
    path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {path}  ({len(rows)} rows)")

# ── Generate all four files ───────────────────────────────────────────────────
for diet in ["SY10", "SY20"]:
    for sex in ["F", "M"]:
        tag  = diet[2:]          # "10" or "20"
        rows = make_rows(diet, sex)
        write_design(OUTDIR / f"AGE_SY{tag}_{sex}.test.txt", rows)

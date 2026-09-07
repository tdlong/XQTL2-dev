#!/bin/bash
# make_patch.sh -- just the passages that changed, as .docx, for pasting into
# the Longevity document a paragraph at a time.
#
#   bash temp_aging/make_patch.sh
#
# Upload to Drive, open with Google Docs, replace the matching paragraphs.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=temp_aging/Longevity_methods_patch.docx
MD=$(mktemp -t patch).md
python3 - "$MD" <<'PY'
import re, sys, subprocess, datetime
m = open('temp_aging/METHODS.md').read()
W = open(sys.argv[1], 'w').write
sha = subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip()
W("# Methods — three passages to replace\n\n")
W("From XQTL2-dev commit %s, %s. Each block replaces the matching passage in\n"
  "\"Longevity\". Nothing else in the Methods changed.\n\n"
  % (sha, datetime.date.today().strftime('%-d %b %Y')))

def para(start):
    i = m.index(start)
    j = m.index("\n\n", i)
    return m[i:j]

W("## 1. Wald test — one sentence\n\n")
W("Replaces the sentence beginning \"on seven degrees of freedom\". "
  "Defines the symbol used below.\n\n")
W(m[m.index("on $\\nu = 7$ degrees"):m.index("Significance is reported")].strip() + "\n\n")

W("## 2. The relationship between the Wald test and h²\n\n")
W("Replaces the whole paragraph. Two corrections: the GWAS analogy carries its\n"
  "own degrees of freedom, and the X-QTL relation carries nu. The old form said\n"
  "E[T] was the non-centrality; for a non-central chi-square E[T] = nu + lambda.\n\n")
W(para("**The relationship between the Wald test") + "\n\n")

W("## 3. Power — first sentence\n\n")
W("Replaces \"E[T] is the non-centrality parameter for the X-QTL design.\", which\n"
  "was not true.\n\n")
W(para("The non-centrality is therefore") + "\n")
PY
pandoc "$MD" --reference-doc=temp_aging/reference.docx -o "$OUT"
python3 temp_aging/fix_math.py "$OUT"
rm -f "$MD"
echo "Wrote $OUT"

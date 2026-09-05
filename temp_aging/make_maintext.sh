#!/bin/bash
# make_maintext.sh -- Results and computational Methods as .docx, for moving
# into the Longevity document.
#
#   bash temp_aging/make_maintext.sh
#
# Upload to Drive, open with Google Docs, copy the blocks across. Equations
# arrive native and editable.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=temp_aging/Longevity_maintext_updates.docx
MD=$(mktemp -t maintext).md
{
  echo "# TRANSFER — for pasting into \"Longevity\""
  echo
  echo "Generated $(date '+%-d %b %Y') from XQTL2-dev, commit $(git rev-parse --short HEAD)."
  echo
  echo "# BLOCK 1 — Results"
  echo
  echo "**Replaces the four Results paragraphs currently in \"Longevity\".**"
  echo
  cat temp_aging/RESULTS.md
  echo
  echo "# BLOCK 2 — Methods, computational sections"
  echo
  echo "**Additive.** The Methods section of \"Longevity\" is a wet-lab outline and"
  echo "contains none of this."
  echo
  sed '1s/^# Methods (incomplete)$/## Methods/' temp_aging/METHODS.md
} > "$MD"
pandoc "$MD" --reference-doc=temp_aging/reference.docx -o "$OUT"
python3 temp_aging/fix_math.py "$OUT"
rm -f "$MD"
echo "Wrote $OUT"

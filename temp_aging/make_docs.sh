#!/bin/bash
# make_docs.sh -- the two documents that go to Drive.
#
#   bash temp_aging/make_docs.sh
#
#   temp_aging/Methods.docx      the computational Methods, whole
#   temp_aging/Supplement.docx   Table S1, the S1/S2 legends, Supplementary Result X
#
# Upload both to Drive, open with Google Docs, copy across.
set -euo pipefail
cd "$(dirname "$0")/.."

sed '1s/^# Methods (incomplete)$/# Methods/' temp_aging/METHODS.md > /tmp/m.md
pandoc /tmp/m.md --reference-doc=temp_aging/reference.docx -o temp_aging/Methods.docx
python3 temp_aging/fix_math.py temp_aging/Methods.docx
echo "Wrote temp_aging/Methods.docx"

bash temp_aging/make_supplement.sh >/dev/null
mv -f temp_aging/Longevity_supplement.docx temp_aging/Supplement.docx
echo "Wrote temp_aging/Supplement.docx"
rm -f /tmp/m.md

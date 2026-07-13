#!/usr/bin/env bash
# diagnose_cluster_state.sh — read-only. Figures out exactly what the cluster's
# dirty working tree is, whether it's real work or just an un-pulled re-sort,
# and whether the xlsx / generator chain is intact. Writes NOTHING to the repo.
#
# Run on the cluster from the repo root:
#   bash scripts_oneoffs/AGE_SY/common/diagnose_cluster_state.sh 2>&1 | tee ~/diag.txt
# then paste ~/diag.txt back.

set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || { echo "not in a git repo"; exit 1; }

line(){ printf '\n========== %s ==========\n' "$1"; }

line "WHERE / BRANCH"
pwd
echo "HEAD:        $(git rev-parse --short HEAD)"
echo "HEAD subject: $(git log -1 --pretty=%s)"

line "FETCH + HOW FAR BEHIND origin/main"
git fetch origin main 2>&1 | sed 's/^/  /' || echo "  (fetch failed)"
echo "origin/main: $(git rev-parse --short origin/main 2>/dev/null || echo '?')"
behind=$(git rev-list --count HEAD..origin/main 2>/dev/null || echo '?')
ahead=$(git rev-list --count origin/main..HEAD 2>/dev/null || echo '?')
echo "commits behind origin/main: $behind"
echo "commits ahead  of origin/main: $ahead"
echo "--- commits the cluster is MISSING (HEAD..origin/main): ---"
git log --oneline HEAD..origin/main 2>/dev/null | sed 's/^/  /' || echo "  (none/err)"

line "GIT STATUS (short)"
git status --short

# --- per-file verdict helper -------------------------------------------------
# Compares a modified tracked file against HEAD, ignoring row order.
# Empty order-insensitive diff => file only re-ordered (safe to discard/take committed).
verdict_reorder_only(){
  local f="$1"
  if ! git ls-files --error-unmatch "$f" >/dev/null 2>&1; then
    echo "  [$f] not tracked — skip"; return; fi
  if git diff --quiet -- "$f"; then
    echo "  [$f] UNMODIFIED"; return; fi
  local content_diff
  content_diff=$(diff <(git show "HEAD:$f" 2>/dev/null | sort) <(sort "$f"))
  local ns; ns=$(git diff --numstat -- "$f" | awk '{print "+"$1" -"$2}')
  if [ -z "$content_diff" ]; then
    echo "  [$f] REORDER-ONLY (content identical to HEAD, just row order)  numstat: $ns  -> SAFE to take committed"
  else
    echo "  [$f] REAL CONTENT DIFFERENCE  numstat: $ns  -> DO NOT blind-discard"
    echo "$content_diff" | sed 's/^/        /' | head -40
    local n; n=$(echo "$content_diff" | wc -l)
    [ "$n" -gt 40 ] && echo "        ... ($n content-diff lines total, truncated)"
  fi
}

line "PER-FILE VERDICT: is each modified file real work, or just an un-pulled re-sort?"
for f in \
  helpfiles/AGE_SY/AGE_SY.bams \
  helpfiles/AGE_SY/AGE_SY10_F.test.txt \
  helpfiles/AGE_SY/AGE_SY10_M.test.txt \
  helpfiles/AGE_SY/AGE_SY20_F.test.txt \
  helpfiles/AGE_SY/AGE_SY20_M.test.txt ; do
  verdict_reorder_only "$f"
done

line "XLSX INTEGRITY (data source for design files)"
if git diff --quiet -- helpfiles/AGE_SY/summary_info_v1.xlsx; then
  echo "  summary_info_v1.xlsx is UNMODIFIED vs HEAD  (cluster copy == committed copy)  OK"
else
  echo "  !! summary_info_v1.xlsx IS MODIFIED vs HEAD — data source diverged, inspect before anything"
fi
echo "  working-copy sha1:  $(git hash-object helpfiles/AGE_SY/summary_info_v1.xlsx)"
echo "  HEAD       sha1:  $(git rev-parse HEAD:helpfiles/AGE_SY/summary_info_v1.xlsx)"
echo "  origin/main sha1: $(git rev-parse origin/main:helpfiles/AGE_SY/summary_info_v1.xlsx 2>/dev/null || echo '?')"

line "HAPLOTYPE PARAMETERS diff (full — it's tiny)"
git --no-pager diff -- helpfiles/AGE_SY/AGE_SY_haplotype_parameters.R | sed 's/^/  /'

line "GENERATOR SANITY (does the committed generator read the xlsx and see R12?)"
echo "  generator input line:"
grep -nE 'XLSX *=' scripts_oneoffs/AGE_SY/common/make_AGE_SY_design_files.py | sed 's/^/    /'
echo "  R12 rows visible to python in the cluster xlsx:"
python3 - <<'PY' 2>&1 | sed 's/^/    /'
try:
    import openpyxl
except Exception as e:
    print("openpyxl not importable here:", e); raise SystemExit
wb = openpyxl.load_workbook("helpfiles/AGE_SY/summary_info_v1.xlsx", data_only=True)
for t in ("Control_Flies","Aged_Flies"):
    ws = wb[t]
    for r in ws.iter_rows(values_only=True):
        rep = r[0]
        if isinstance(rep,str) and rep.strip() in ("Rep12","R12","12"):
            rr=[x for x in r if x is not None]
            print(t, rr)
PY

line "NOTABLE UNTRACKED (informational only — not touched)"
git status --short | grep '^??' | sed 's/^/  /' || true

line "DONE"
echo "Report ready. Paste this whole output back."

#!/usr/bin/env bash
# test_3rd_allele_cap.chrX.sh — test the "rare 3rd allele + -d 1000 cap" idea.
#
# Hypothesis (chrX only, ~5x faster than whole genome):
#   The SNPs the two callers disagree on are sites where
#     (1) true depth > 1000 in some sample (so -d 1000 downsamples),
#     (2) a real 3rd allele exists (site is multiallelic when uncapped),
#     (3) capping to 1000 flips the call between biallelic (kept by -m2 -M2)
#         and multiallelic (dropped) — which is why -t and -r, keeping different
#         1000-read subsets, disagree.
#
# For the chrX differing sites AND a control set of agreeing sites, this runs
# bcftools mpileup|call twice — uncapped (-d 100000) and capped (-d 1000) — and
# reports, per site: number of ALT alleles each way, and the true max depth.
#
# Run from repo root on the cluster:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/test_3rd_allele_cap.chrX.sh
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true

CHR=chrX
REF=pipeline/ref/dm6.fa
BAMS=helpfiles/AGE_SY/AGE_SY.bams
a=process/AGE_SY_v3/RefAlt.$CHR.txt
b=process/AGE_SY_v3_tiled/RefAlt.$CHR.txt
for f in "$a" "$b" "$REF" "$BAMS"; do [[ -e "$f" ]] || { echo "missing: $f"; exit 1; }; done

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# positions where the two callers disagree (present in only one)
awk 'NR==FNR{ if(FNR>1)A[$2]=1; next } FNR>1{ B[$2]=1; if(!($2 in A)) print "'"$CHR"'\t"$2 }
     END{ for(p in A) if(!(p in B)) print "'"$CHR"'\t"p }' "$a" "$b" | sort -k2,2n > "$TMP/diff.pos"
# control: agreeing positions (in both), random 200
awk 'NR==FNR{ if(FNR>1)A[$2]=$0; next } FNR>1 && ($2 in A) && A[$2]==$0 { print "'"$CHR"'\t"$2 }' "$a" "$b" \
  | shuf | head -200 | sort -k2,2n > "$TMP/ctrl.pos"

# per-site allele count + max depth, at a given -d, for a position list
# emits: pos  nALT  maxDP
call_at () {  # $1 = -d value   $2 = positions file
  bcftools mpileup -I -d "$1" -a FORMAT/DP -f "$REF" -b "$BAMS" -R "$2" 2>/dev/null \
   | bcftools call -m 2>/dev/null \
   | bcftools query -f '%POS\t%ALT[\t%DP]\n' \
   | awk '{ nalt=($2=="."?0:split($2,aa,",")); mx=0; for(i=3;i<=NF;i++){d=($i=="."?0:$i+0); if(d>mx)mx=d}
            print $1"\t"nalt"\t"mx }'
}

analyze () {  # $1 = positions file   $2 = label
  call_at 100000 "$1" | sort -k1,1n > "$TMP/uncap"
  call_at 1000   "$1" | sort -k1,1n > "$TMP/cap"
  join -1 1 -2 1 "$TMP/uncap" "$TMP/cap" \
    | awk -v lab="$2" '
        # fields: pos  naltU maxDP  naltC maxDPc
        { n++
          if($3>1000) deep++
          if($2>=2) multiU++
          if($2!=$4) flip++
          if($3>1000 && $2>=2) both++ }
        END{ if(n==0){print "  "lab": (no sites)"; next}
          printf "  %-22s %3d sites | depth>1000: %3d (%.0f%%) | multiallelic-uncapped: %3d (%.0f%%) | cap-changes-#alleles: %3d (%.0f%%) | depth>1000 AND 3rd-allele: %3d (%.0f%%)\n",
                 lab, n, deep,100*deep/n, multiU,100*multiU/n, flip,100*flip/n, both,100*both/n }'
}

echo "chrX: differing sites = $(wc -l < "$TMP/diff.pos"),  control sites = $(wc -l < "$TMP/ctrl.pos")"
echo "(nALT: number of ALT alleles bcftools calls; >=2 means a 3rd allele => dropped by -m2 -M2)"
echo
analyze "$TMP/diff.pos" "DIFFERING sites:"
analyze "$TMP/ctrl.pos" "CONTROL sites:"
echo
echo "Example differing sites (pos | nALT_uncapped | trueMaxDepth | nALT_capped):"
call_at 100000 "$TMP/diff.pos" | sort -k1,1n > "$TMP/u2"
call_at 1000   "$TMP/diff.pos" | sort -k1,1n > "$TMP/c2"
join "$TMP/u2" "$TMP/c2" | awk '{printf "  %-10s nALTuncap=%s  trueMaxDP=%-7s nALTcap=%s\n",$1,$2,$3,$5}' | head -25

#!/bin/bash
#SBATCH --job-name=vet_founders
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/vet_founders_%j.out
###############################################################################
# vet_founders.sh — STAGE 1 of the one-SNP dissection. For each candidate
# (v3-vs-v4 disagreement) SNP, look it up in the NEW annotated founder catalog
# (catalog.annot.tsv.gz: CHROM POS REF ALT dist_indel AD_<founder>=ref,alt) and
# decide which SNPs the FOUNDERS say are clean: every founder fixed (alt freq
# <=0.03 or >=0.97), both REF and ALT represented, and far from a founder indel.
# Those are the sites where the "true" answer is unambiguous, so a v3/v4 count
# disagreement there is purely a read-handling artifact worth dissecting.
#
# Also cuts a small BAM slice per (sample,pos) for Stage 2 (local alignment view).
#
#   sbatch scripts_oneoffs/AGE_SY/catalog_v5/vet_founders.sh
###############################################################################
set -uo pipefail
module load bcftools/1.21 2>/dev/null || true
module load samtools/1.21 2>/dev/null || true
ANNOT=process/AGE_SY_v5/Catalog/catalog.annot.tsv.gz
CAND=scripts_oneoffs/AGE_SY/catalog_v5/dissect_candidates.tsv
SLICEDIR=process/AGE_SY_v5/slices
BAMDIR=data/bam/AGE_SY
for f in "$ANNOT" "$CAND"; do [[ -s "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done
mkdir -p "$SLICEDIR"

echo "############ STAGE 1: founder vet of $(($(wc -l < "$CAND")-1)) candidate SNPs ############"
echo "founder AD-column order: $(zcat "$ANNOT" | head -1 | cut -f6- | sed 's/AD_//g')"
echo "CLEAN = all founders fixed (<=0.03 or >=0.97), >=1 REF & >=1 ALT, dist_indel>=25"
echo

awk -F'\t' -v OFS='\t' '
  NR==FNR { if(FNR>1) want[$2]=1; next }            # candidate positions
  FNR==1  { for(i=6;i<=NF;i++){ n=$i; sub(/^AD_/,"",n); fn[i]=n } next }   # annot header -> founder names
  ($2 in want) {
    seen[$2]=1; d=$5; nref=0; nalt=0; nint=0; na=0; vec="";
    for(i=6;i<=NF;i++){
      split($i,a,","); tot=a[1]+a[2];
      if(tot==0){ na++; vec=vec sprintf("%s=NA ",fn[i]); continue }
      f=a[2]/tot;
      if(f<=0.03) nref++; else if(f>=0.97) nalt++; else nint++;
      vec=vec sprintf("%s=%.2f ",fn[i],f);
    }
    clean=(nint==0 && nref>=1 && nalt>=1 && d>=25) ? "CLEAN" : "---  ";
    printf "%s  %s:%d  dist_indel=%-6s  Rfix=%d Afix=%d inter=%d na=%d   %s\n",
           clean,$1,$2,d,nref,nalt,nint,na,vec;
  }
  END { }
' "$CAND" <(zcat "$ANNOT") | sort

echo
echo "--- candidates NOT found in v5 annotated catalog (dropped as non-biallelic in v5?) ---"
awk -F'\t' 'NR==FNR{if(FNR>1)c[$2]=$1; next} FNR>1{seen[$2]=1} END{for(p in c) if(!(p in seen)) print c[p]":"p}' \
  "$CAND" <(zcat "$ANNOT") | sort

echo
echo "############ cutting BAM slices (chr:pos +/-300) for Stage 2 ############"
tail -n +2 "$CAND" | while IFS=$'\t' read -r chr pos samp rest; do
  bam="$BAMDIR/$samp.bam"
  [[ -s "$bam" ]] || { echo "  no bam: $bam (skip)"; continue; }
  out="$SLICEDIR/${samp}.${chr}_${pos}.bam"
  samtools view -b "$bam" "$chr:$((pos-300))-$((pos+300))" > "$out" 2>/dev/null && samtools index "$out" \
    && echo "  $out  ($(samtools view -c "$out") reads)"
done
echo "slices -> $SLICEDIR  (rsync local for interactive tview)"
echo "############ END STAGE 1 ############"

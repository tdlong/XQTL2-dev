#!/bin/bash
#SBATCH --job-name=snp_loss_grid
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:20:00
#SBATCH -o logs/AGE_SY/snp_loss_grid.out
# STEP 2: SNP-loss table for the 12 filter combos (indel 0/10/20/50 x maxaf 3/5/7%),
# holding >=10 reads + polymorphic. Pure count off catalog.annot.tsv.gz (BAQ-on build).
# Also emits <catdir>/snp_pass.tsv.gz (CHROM POS dist_indel pass03 pass05 pass07) for the
# concordance step. B5 ignored on chr2L.   ->  logs/AGE_SY/snp_loss_grid.out
#   sbatch scripts_oneoffs/AGE_SY/catalog_v6/snp_loss_grid.sh <catdir>
set -uo pipefail
CATDIR=${1:-process/AGE_SY_v6/Catalog}
ANNOT="$CATDIR/catalog.annot.tsv.gz"
PASS="$CATDIR/snp_pass.tsv.gz"
[[ -s "$ANNOT" ]] || { echo "missing $ANNOT" >&2; exit 1; }

zcat "$ANNOT" | awk -F'\t' -v PASS="gzip > $PASS" '
BEGIN{
  ng=split("0 10 20 50",G," ");
  nm=split("0.03 0.05 0.07",M," ");
  print "CHROM\tPOS\tdist_indel\tpass03\tpass05\tpass07" | PASS;
}
NR==1{ for(i=6;i<=NF;i++){n=$i; sub(/^AD_/,"",n); if(n=="B5") b5=i} next }
{
  chr=$1; d=$5+0; cand++;
  depthok=1;
  for(k=1;k<=nm;k++){ nf[k]=1; hr[k]=0; ha[k]=0 }
  for(i=6;i<=NF;i++){
    if(chr=="chr2L" && i==b5) continue;                # B5 exempt on chr2L
    split($i,a,","); dep=a[1]+a[2];
    if(dep<10) depthok=0;
    af=(dep>0? a[2]/dep : -1);
    for(k=1;k<=nm;k++){ m=M[k];
      if(af<=m) hr[k]=1; else if(af>=1-m) ha[k]=1; else nf[k]=0 }
  }
  p3=(depthok && nf[1] && hr[1] && ha[1])?1:0;
  p5=(depthok && nf[2] && hr[2] && ha[2])?1:0;
  p7=(depthok && nf[3] && hr[3] && ha[3])?1:0;
  print chr"\t"$2"\t"d"\t"p3"\t"p5"\t"p7 | PASS;
  # loss-table tallies: pass founder rules @m AND dist_indel>=g
  pk[1]=p3; pk[2]=p5; pk[3]=p7;
  for(k=1;k<=nm;k++){ if(!pk[k]) continue; for(g=1;g<=ng;g++) if(d>=G[g]) cnt[g,k]++ }
}
END{
  close(PASS);
  printf "candidate SNPs (biallelic in founders, BAQ-on build): %d\n\n", cand;
  printf "SNPs kept (rows = >= bp from indel, cols = maxaf; >=10 reads + polymorphic held)\n";
  printf "%-10s", "indel>=";
  for(k=1;k<=nm;k++) printf "%12s", "af"M[k];
  print "";
  for(g=1;g<=ng;g++){ printf "%-10s", G[g];
    for(k=1;k<=nm;k++) printf "%12d", cnt[g,k];
    print "" }
  printf "\nper-SNP pass table -> %s\n", "'"$PASS"'";
}'

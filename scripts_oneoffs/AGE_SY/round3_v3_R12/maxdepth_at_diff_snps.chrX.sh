#!/usr/bin/env bash
# maxdepth_at_diff_snps.chrX.sh
#
# The guess, made from data we already have (no mpileup, seconds of awk):
# the two callers disagree only on WHICH SNPs pass (shared SNPs have identical
# counts). For each DISAGREEING SNP, look PER SAMPLE and ask: is one sample
# sitting near the -d 1000 ceiling (i.e. that sample was subsampled, so the
# chunked vs whole stream could have kept different reads -> different QUAL)?
#
# Both tables are built at -d 1000, so every sample's ref+alt is ceilinged at
# ~1000 (post quality-filter, so a capped sample reads ~700-950, not exactly
# 1000). We therefore show the DISTRIBUTION of max per-sample depth, not a hard
# threshold:
#   - differing SNPs pile up near the ceiling  => the cap drives the disagreement
#     (consistent with subsampling).
#   - differing SNPs with max per-sample depth well BELOW the ceiling (say <500,
#     raw < ~1000 for any plausible hq fraction) => no sample capped => both
#     callers saw identical reads => the disagreement is NOT subsampling.
# Compare against the AGREEING (shared) SNPs as a baseline.
#
# Reads only; run on the login node:
#   bash scripts_oneoffs/AGE_SY/round3_v3_R12/maxdepth_at_diff_snps.chrX.sh
set -uo pipefail

A=process/AGE_SY_v3/RefAlt.chrX.txt          # OLD = whole-chromosome caller (-t)
B=process/AGE_SY_v3_tiled/RefAlt.chrX.txt     # NEW = tiled caller (-r, 5 Mb chunks)
for f in "$A" "$B"; do [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }; done

awk '
  # max per-sample depth (ref+alt) for one body line -> "max<TAB>argmax_k<TAB>n>=900"
  function summarize(line,   nf,f,k,ns,d,mx,mi,nn){
    nf=split(line,f," "); ns=(nf-2)/2; mx=0; mi=0; nn=0
    for(k=1;k<=ns;k++){ d=f[2*k+1]+f[2*k+2]; if(d>mx){mx=d;mi=k}; if(d>=900)nn++ }
    return mx"\t"mi"\t"nn
  }
  function bucket(m){ return (m<300)?"<300":(m<500)?"300-499":(m<700)?"500-699":(m<900)?"700-899":">=900" }

  # ---- file A: header (sample names) then body ----
  FNR==1 && NR==FNR { n=split($0,h," "); for(k=1;k<=(n-2)/2;k++){ nm=h[2*k+1]; sub(/^REF_/,"",nm); name[k]=nm } next }
  NR==FNR { s=summarize($0); split(s,z,"\t"); Amx[$2]=z[1]; Asm[$2]=z[2]; Ann[$2]=z[3]; next }
  # ---- file B: skip header, then body ----
  FNR==1 { next }
  { s=summarize($0); split(s,z,"\t"); Bmx[$2]=z[1]; Bsm[$2]=z[2]; Bnn[$2]=z[3] }

  END{
    # classify
    for(p in Amx){ if(p in Bmx){ shared[p]=1 } else { oldonly[p]=1 } }
    for(p in Bmx){ if(!(p in Amx)) newonly[p]=1 }
    ns=0; for(p in shared)ns++; no=0; for(p in oldonly)no++; nn=0; for(p in newonly)nn++

    print  "==================================================================="
    print  "MAX PER-SAMPLE DEPTH at disagreeing SNPs  (both tables -d 1000)"
    print  "  A/OLD = whole-chr (-t)   B/NEW = tiled (-r)"
    print  "==================================================================="
    printf "shared: %d    OLD-only: %d    NEW-only: %d    (disagree: %d)\n\n", ns,no,nn,no+nn

    # distribution: disagreeing vs shared
    print  "distribution of MAX per-sample depth (what fraction has a sample near the ceiling):"
    printf "  %-10s %8s %8s %8s %8s %8s\n","group","<300","300-499","500-699","700-899",">=900"
    delta=0
    # disagreeing
    for(p in oldonly){ b[bucket(Amx[p])]++; d0++ }
    for(p in newonly){ b[bucket(Bmx[p])]++; d0++ }
    printf "  %-10s %8d %8d %8d %8d %8d\n","disagree", b["<300"]+0,b["300-499"]+0,b["500-699"]+0,b["700-899"]+0,b[">=900"]+0
    # shared baseline
    delete b
    for(p in shared){ b[bucket(Amx[p])]++ }
    printf "  %-10s %8d %8d %8d %8d %8d\n","shared", b["<300"]+0,b["300-499"]+0,b["500-699"]+0,b["700-899"]+0,b[">=900"]+0

    # conservative residual: disagreeing SNPs with NO sample plausibly capped (max<500)
    print  ""
    print  "CONSERVATIVE RESIDUAL — disagreeing SNPs with max per-sample depth < 500"
    print  "  (raw < ~1000 for any plausible hq fraction => no sample capped => NOT subsampling):"
    printf "  %-12s %-8s %-8s %-20s\n","pos","caller","maxDP","sample_at_max"
    res=0
    for(p in oldonly){ if(Amx[p]<500){ printf "  chrX:%-8s %-8s %-8d %-20s\n",p,"OLD",Amx[p],name[Asm[p]]; res++ } }
    for(p in newonly){ if(Bmx[p]<500){ printf "  chrX:%-8s %-8s %-8d %-20s\n",p,"NEW",Bmx[p],name[Bsm[p]]; res++ } }
    printf "  --> %d of %d disagreeing SNPs have NO plausibly-capped sample.\n", res, no+nn
    if(res>0) print  "      => those cannot be a -d subsampling artifact. Not everything is sampling."
    if(res==0)print  "      => every disagreeing SNP has a sample near the ceiling. Consistent with subsampling."

    # a few high-maxDP examples, to eyeball the near-ceiling cases
    print  ""
    print  "sample of near-ceiling disagreeing SNPs (max>=900), to eyeball:"
    c=0
    for(p in oldonly){ if(Amx[p]>=900 && c<10){ printf "  chrX:%-8s OLD  maxDP=%-6d  %s\n",p,Amx[p],name[Asm[p]]; c++ } }
    for(p in newonly){ if(Bmx[p]>=900 && c<20){ printf "  chrX:%-8s NEW  maxDP=%-6d  %s\n",p,Bmx[p],name[Bsm[p]]; c++ } }
  }
' "$A" "$B"

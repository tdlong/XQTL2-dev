review for XQTL2-dev

Development version of XQTL2 pipeline. HPC has raw and intermediates.
Per your correction, process/ is HPC-only (not cached on desktop).

Per-folder classifications:

scripts: 111
scripts_oneoffs: 111
helpfiles: 111
logs: 111

process: 100

Extra .gitignore lines:

# specific loose files to exclude/delete
R.haps.chr2L.out.rds
repro_error.tar.gz

# SLURM and pipeline output patterns (from prior .gitignore, preserved)
slurm*
*.out
*_report.txt

# big binary data
*.bam
*.bai
*.rds
*.RDS
*.tar
*.tar.gz

# defensive: directories that would be big if created
data/
ref/
figures/
output/
configs/
process.*/

# loose images at project root (real ones live in figures/)
*.png
*.pdf
*.jpg
*.jpeg

# machine-specific pipeline symlink
pipeline

# specific known file
check_cluster_state.out.txt

# override: logs/ IS tracked (your workflow - Claude reads them for debugging)
!logs/
!logs/**

Notes:

logs/ stays in git so Claude can debug from any machine.

process/ is HPC-canonical only. Its desktop copy is currently empty
(any content there was created locally or rsynced from HPC, never git).

R.haps.chr2L.out.rds and repro_error.tar.gz are excluded from git
because you asked them deleted.

HPC path: /dfs7/adl/tdlong/fly_pool/XQTL2/

Files to delete on your Mac (device_bash cannot rm):

  rm /Volumes/Raw_store/projects/XQTL2-dev/R.haps.chr2L.out.rds
  rm /Volumes/Raw_store/projects/XQTL2-dev/repro_error.tar.gz

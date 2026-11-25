# Haplotype parameters for the full 2022 malathion experiment: 8 pools, pA founders.
# Same window settings as helpfiles/malathion/hap_params.R; only names_in_bam differs.
# Rep-2 names follow the rep-1 pattern and are checked against the RefAlt header
# by scripts_oneoffs/malathion/rerun_2022.sh before anything is submitted.
founders <- c("A1","A2","A3","A4","A5","A6","A7","AB8")
names_in_bam <- c("c_F_1","c_M_1","s_F_1","s_M_1","c_F_2","c_M_2","s_F_2","s_M_2")
step <- 5000
size <- 50000
h_cutoff <- 2.5

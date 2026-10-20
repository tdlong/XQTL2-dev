# palette.R -- colours shared across the AGE_SY figures.  source() this; do not
# copy the values, or "the same colour scheme" becomes a coincidence that a later
# edit quietly breaks.
#
#   source("temp_aging/analysis/palette.R")
#
# Founder colours are the XQTL2.Xplore palette (github.com/tdlong/XQTL2.Xplore),
# assigned by sorted founder name so a founder is the same colour in every panel
# of every figure.  make_figure1.R and make_figure2.R still define these inline;
# the values are identical and those two should be migrated when next touched.

FOUNDER_COL <- c(AB8 = "#E69F00", B1 = "#56B4E9", B2 = "#009E73", B3 = "#F0E442",
                 B4  = "#D55E00", B5 = "#CC79A7", B6 = "#000000", B7 = "#990099")

TRT_COL <- c("SY10 female" = "#F49AC2", "SY20 female" = "#D62728",
             "SY10 male"   = "#8EC7E8", "SY20 male"   = "#1F4E9C")

TRT_LAB <- c(SY10_F = "SY10 female", SY20_F = "SY20 female",
             SY10_M = "SY10 male",   SY20_M = "SY20 male")

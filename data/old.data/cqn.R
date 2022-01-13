library(cqn)
library(tidyverse)

peaklist <- read.delim(file = "/public/groups/forsberglab/atac_analysis/200519_bamscale/final_peaks_0622.txt")
counttable <- read.delim(file = "/public/groups/forsberglab/atac_analysis/200519_bamscale/allcounts_0622.txt", row.names = 1)

atac_cqn_peaklist0622 <- cqn(counttable, peaklist$bias, peaklist$width, sizeFactors = NULL, lengthMethod = "fixed")

save(atac_cqn_peaklist0622, file = "/public/groups/forsberglab/atac_analysis/200519_bamscale/cqn_peaklist_0622.Rdata")

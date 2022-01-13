library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
library(BSgenome.Mmusculus.UCSC.mm10)
register(SerialParam())

load("/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/filtered_counts_211129.Rdata")

peaks_dev_211129 <- computeDeviations(object = peak_counts_filtered_211129)

save(peaks_dev_211129, file = "/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/peaks_dev_211129.Rdata")

library(ggplot2)

peaks_variability_211129 <- computeVariability(peaks_dev_211129)

save(peaks_variability_211129, file = "/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/variability_211129.Rdata")

rm(peaks_variability_211129)

tsne_peaks_211129 <- deviationsTsne(peaks_dev_211129, threshold = 1.5, perplexity = 10, shiny = FALSE)

save(tsne_peaks_211129, file = "/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/tsne_211129.Rdata")

rm(tsne_peaks_211129)

cor_211129 <- getSampleCorrelation(peaks_dev_211129)

save(cor_211129, file = "/public/groups/forsberglab/211129_chromvar_pros/cor_211129.Rdata")


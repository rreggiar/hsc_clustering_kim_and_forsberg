library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
register(MulticoreParam(8))
library(BSgenome.Mmusculus.UCSC.mm10)

peakfile <-"/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/idr_hsc_pros_211129_peaks_sorted.bed"

peaks <- readNarrowpeaks(peakfile, width = 500, non_overlapping = FALSE)
peaks

clp1 <- "/public/groups/forsberglab/atac_libraries/clp/ATAC34_CLP_CFEM026_S9_R1_001.merged.nodup.bam"
clp2 <- "/public/groups/forsberglab/atac_libraries/clp/cfjk017_R1.trim.merged.nodup.bam"
cmp1 <- "/public/groups/forsberglab/atac_libraries/cmp/em013_cat_R1.trim.merged.nodup.bam"
cmp2 <- "/public/groups/forsberglab/atac_libraries/cmp/jk004_cat_R1.trim.merged.nodup.bam"
gmp1 <- "/public/groups/forsberglab/atac_libraries/gmp/em012_cat_R1.trim.merged.nodup.bam"
gmp2 <- "/public/groups/forsberglab/atac_libraries/gmp/ATAC33_GMP_CFEM025_S8_R1_001.merged.nodup.bam"
hsc1 <- "/public/groups/forsberglab/atac_libraries/hsc/jk009_allcat_R1.merged.nodup.bam"
hsc2 <- "/public/groups/forsberglab/atac_libraries/hsc/rs003.merged.nodup.bam"
mep1 <- "/public/groups/forsberglab/atac_libraries/mep/em014_cat_R1.merged.nodup.bam"
mep2 <- "/public/groups/forsberglab/atac_libraries/mep/em020.merged.nodup.bam"
mpp1 <- "/public/groups/forsberglab/atac_libraries/mpp/ATAC33_MPP_CFEM024_S7_R1_001.merged.nodup.bam"
mpp2 <- "/public/groups/forsberglab/atac_libraries/mpp/rs004.merged.nodup.bam"
prob1 <- "/public/groups/forsberglab/atac_libraries/prob/ATAC35_ProB_CFEM029_S10_R1_001.trim.merged.nodup.no_chrM_MT.bam"
prob2 <- "/public/groups/forsberglab/atac_libraries/prob/jk006_cat_R1.trim.merged.nodup.no_chrM_MT.bam"
prot1 <- "/public/groups/forsberglab/atac_libraries/prot/em007_cat_R1.trim.merged.nodup.no_chrM_MT.bam"
prot2 <- "/public/groups/forsberglab/atac_libraries/prot/jk007_cat_R1.trim.merged.nodup.no_chrM_MT.bam"

inbam <-c(hsc1,hsc2,mpp1,mpp2,cmp1,cmp2,gmp1,gmp2,mep1,mep2,clp1,clp2,prob1,prob2,prot1,prot2)
peak_counts <- getCounts(inbam, peaks, by_rg = FALSE, format = "bam", paired = TRUE,
colData = DataFrame(celltype = c("HSC","HSC","MPP","MPP","CMP","CMP","GMP","GMP","MEP","MEP","CLP","CLP","ProB","ProB","ProT","ProT")))

peak_counts <- addGCBias(peak_counts, genome = BSgenome.Mmusculus.UCSC.mm10)

peak_counts <- sort(peak_counts)

peak_counts_filtered_211129 <- filterPeaks(peak_counts, non_overlapping = TRUE)

write.table(rowRanges(peak_counts_filtered_211129),file="/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/hsc_pros_peaklist_211129.txt",quote= FALSE,sep = "\t")

save(peak_counts_filtered_211129,file = "/public/groups/forsberglab/atac_analysis/211129_chromvar_pros/filtered_counts_211129.Rdata")


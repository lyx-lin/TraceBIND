projectName = "NR3C1/Control6"

print(projectName)

library(MASS)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)
library(tensorflow)
# library(Signac)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

if(!dir.exists(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/", projectName))){
  system(paste("mkdir -p", paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/", projectName)))
}

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/", projectName))

tf = readxl::read_excel("~/nzhanglab/project/linyx/footprints/PRINT/data/NR3C1/Control6/41467_2022_32972_MOESM14_ESM.xlsx",
                        col_names = FALSE)
colnames(tf) = c('chr', 'start', 'end')

tf = tf %>%
  subset(chr %in% paste0('chr', 1:22))
write.table(tf, 
            "~/nzhanglab/project/linyx/footprints/PRINT/data/NR3C1/hTERT_RPTEC_GR_CUT&RUN_filtering.txt")

peaks = read.table('~/nzhanglab/data/ParkerWilson/kidney_multiome_atlas/peaks/peaks.csv', 
                   sep = ',', header = TRUE)

peaks = Signac::StringToGRanges(peaks$Peaks, sep = c(":", "-"))

tf_chip = GenomicRanges::GRanges(seqnames = tf$chr, 
                                 ranges = IRanges::IRanges(start = tf$start, 
                                                           end = tf$end))
pairs = IRanges::findOverlapPairs(tf_chip, peaks)
tf_interest = pairs@first
regions = pairs@second
regionsBed = data.frame("chr" = regions@seqnames, 
                        "start" = as.data.frame(regions@ranges)$start, 
                        "end" = as.data.frame(regions@ranges)$end)
regions = GenomicRanges::GRanges(seqnames = regionsBed$chr, 
                                 ranges = IRanges::IRanges(start = regionsBed$start, 
                                                           end = regionsBed$end))
regions = reduce(regions)
regionsBed = data.frame("chr" = regions@seqnames, 
                        "start" = as.data.frame(regions@ranges)$start, 
                        "end" = as.data.frame(regions@ranges)$end)

tf_interest = data.frame("chr" = tf_interest@seqnames, 
                           "start" = as.data.frame(tf_interest@ranges)$start, 
                           "end" = as.data.frame(tf_interest@ranges)$end)

regionsBed = unique(regionsBed)
tf_interest = unique(tf_interest)

write.table(regionsBed, 
            paste0("peaks_region.txt"))
write.table(tf_interest, 
            paste0("tf_region.txt"))

get_bias(regionsBed, 
         referenceGenome = 'hg38', 
         path = paste0(getwd(), '/freezed_finetuned_'),
         code_path = '../../../../', 
         model_use = 'Tn5_NN_model_Control_6_freezed_finetuned.h5'
         )

metadata = read.csv(paste0("~/nzhanglab/data/ParkerWilson/kidney_multiome_atlas/kidney_chasm/Control_6/bcanno.csv"))
barcodeGroups = metadata[c('barcode', 'celltype')]
frags_path = "~/nzhanglab/data/ParkerWilson/kidney_multiome_atlas/kidney_fragment_files/Control_6.fragments.tsv.gz"
counts = get_count(frags_path, 
                   regionsBed, 
                   barcodeGroups = barcodeGroups)
saveRDS(counts, 'counts.rds')

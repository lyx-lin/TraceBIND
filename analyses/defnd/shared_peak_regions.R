library(GenomicRanges)
library(tidyverse)

peak_v22 = read.table("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-V22-KYC-LN-01-DEFND-ATAC/peaks.bed")
peak_u21_defnd = read.table("~/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-DEFND-ATAC/peaks.bed")
peak_u21_multiome = read.table("~/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-Multiome-ATAC/peaks.bed")
colnames(peak_v22) = c('chr', 'start', 'end')
colnames(peak_u21_defnd) = c('chr', 'start', 'end')
colnames(peak_u21_multiome) = c('chr', 'start', 'end')

overlaps = findOverlapPairs(makeGRangesFromDataFrame(peak_u21_defnd), 
                            makeGRangesFromDataFrame(peak_u21_multiome))
intersection = pintersect(overlaps@first, 
                          overlaps@second)
overlaps = findOverlapPairs(intersection, 
                            makeGRangesFromDataFrame(peak_v22))
intersection = pintersect(overlaps@first, 
                          overlaps@second)

intersection = as.data.frame(intersection)
intersection = intersection %>%
  dplyr::select(seqnames, start, end)
colnames(intersection) = c('chr', 'start', 'end')
intersection = intersection[intersection$end - intersection$start > 200, ]
write.table(intersection, 
            '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/peaks_region.bed')

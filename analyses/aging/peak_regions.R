library(GenomicRanges)
library(tidyverse)

peak_82 = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-Multiome-ATAC/peaks.bed")
peak_56 = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-Q17-KYC-LN-01-Multiome-ATAC/peaks.bed")
peak_30 = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-L12-KYC-LN-01-Multiome-ATAC/peaks.bed")
peak_16 = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/RAGE24-C03-KYC-LN-01-Multiome-ATAC/peaks.bed")

colnames(peak_82) = c('chr', 'start', 'end')
colnames(peak_56) = c('chr', 'start', 'end')
colnames(peak_30) = c('chr', 'start', 'end')
colnames(peak_16) = c('chr', 'start', 'end')

### all -----
peak_82 <- makeGRangesFromDataFrame(peak_82)
peak_56 = makeGRangesFromDataFrame(peak_56)
peak_30 <- makeGRangesFromDataFrame(peak_30)
peak_16 <- makeGRangesFromDataFrame(peak_16)

combined.peaks <- reduce(x = c(peak_16, peak_30, peak_56, peak_82))
combined.peaks = combined.peaks[width(combined.peaks) >= 800]
combined.peaks = as.data.frame(combined.peaks)
combined.peaks = combined.peaks %>%
  dplyr::select(seqnames, start, end)
colnames(combined.peaks) = c('chr', 'start', 'end')
write.table(combined.peaks, 
            "~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/all_peaks_region.bed")

### overlapping -----
overlaps_1 = findOverlapPairs(makeGRangesFromDataFrame(peak_82), 
                              makeGRangesFromDataFrame(peak_56))
intersection_1 = pintersect(overlaps_1@first, 
                            overlaps_1@second)
intersection_1 = intersection_1[width(intersection_1) > 800]
overlaps_2 = findOverlapPairs(intersection_1, 
                              makeGRangesFromDataFrame(peak_30))
intersection_2 = pintersect(overlaps_2@first, 
                            overlaps_2@second)
intersection_2 = intersection_2[width(intersection_2) > 800]

overlaps_3 = findOverlapPairs(intersection_2, 
                              makeGRangesFromDataFrame(peak_16))
intersection_3 = pintersect(overlaps_3@first, 
                            overlaps_3@second)
intersection_3 = intersection_3[width(intersection_3) > 800]

intersection_3 = as.data.frame(intersection_3)
intersection_3 = intersection_3 %>%
  dplyr::select(seqnames, start, end)
colnames(intersection_3) = c('chr', 'start', 'end')
write.table(intersection_3, 
            '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/peaks_region.bed')

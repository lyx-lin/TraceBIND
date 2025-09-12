library(Seurat)
library(Signac)
library(tidyverse)

setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")

seuratobject_whole = readRDS("seuratobject_whole_peaks.rds")
seuratobject_whole = seuratobject_whole %>%
  subset(celltype %in% c("PT", "PT-MT"))

mean_peaks = seuratobject_whole@meta.data %>% 
  group_by(dataset) %>% 
  summarise(mean_peak = mean(nFeature_ATAC))
cor.test(as.numeric(unlist(mean_peaks[2])), 
         c(16, 30, 56, 82))

ggplot(data = mean_peaks,
       aes(x = as.numeric(as.factor(dataset)), y = mean_peak)) + 
  geom_boxplot(data = seuratobject_whole@meta.data, 
               aes(x = as.numeric(as.factor(dataset)), 
                   y = nFeature_ATAC, 
                   group = dataset)) + 
  xlab('ages') + 
  ylab("peaks counts") + 
  geom_point(color = 'red') + 
  geom_smooth(method = lm,
              se = FALSE, color = 'red')

ggsave("nfeatures_ages.pdf", 
       height = 5, width = 6)

peaks_c03 = read.table("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-C03-KYC-LN-01-Multiome-ATAC/peaks.bed")
peaks_l12 = read.table("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-L12-KYC-LN-01-Multiome-ATAC/peaks.bed")
peaks_q17 = read.table("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-Q17-KYC-LN-01-Multiome-ATAC/peaks.bed")
peaks_u21 = read.table("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-Multiome-ATAC/peaks.bed")
sum_peaks_c03 = sum(peaks_c03$V3 - peaks_c03$V2 + 1)
sum_peaks_l12 = sum(peaks_l12$V3 - peaks_l12$V2 + 1)
sum_peaks_q17 = sum(peaks_q17$V3 - peaks_q17$V2 + 1)
sum_peaks_u21 = sum(peaks_u21$V3 - peaks_u21$V2 + 1)
sum_peaks_c03 = sum_peaks_c03/2440000000
sum_peaks_l12 = sum_peaks_l12/2440000000
sum_peaks_q17 = sum_peaks_q17/2440000000
sum_peaks_u21 = sum_peaks_u21/2440000000
cor.test(c(sum_peaks_c03, sum_peaks_l12, sum_peaks_q17, sum_peaks_u21), 
         c(16, 30, 56, 82))
ggplot() + 
  geom_point(aes(x = c(16, 30, 56, 82), 
                 y = 100*c(sum_peaks_c03, sum_peaks_l12, sum_peaks_q17, sum_peaks_u21))) + 
  geom_smooth(aes(x = c(16, 30, 56, 82), 
                  y = 100*c(sum_peaks_c03, sum_peaks_l12, sum_peaks_q17, sum_peaks_u21)), 
              se = F, method = 'lm', color = 'red') + 
  xlab('age') + 
  ylab("% of the genome in accessible peaks")
ggsave("percentage_genome_ages.pdf", 
       height = 5, width = 6)

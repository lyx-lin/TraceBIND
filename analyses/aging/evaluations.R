library(tidyverse)
library(GenomicRanges)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(GGally)
library(hdrcde)
library(ggdensity)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggpubr)
library(ggrepel)
library(gridExtra)
theme_Publication <- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus"))
source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')
### coverage ----
counts_16 = readRDS("age_16_weeks/C03/kidney_multiome/counts.rds")
counts_30 = readRDS("age_30_weeks/L12/kidney_multiome/counts_subset.rds")
counts_56 = readRDS("age_56_weeks/Q17/kidney_multiome/counts_subset.rds")
counts_82 = readRDS("age_82_weeks/U21/kidney_multiome/counts_subset.rds")

average_16 = NULL
for (i in 1:length(counts_16)){
  average_16 = c(average_16, 
                 sum(counts_16[[i]]$count))
}

average_30 = NULL
for (i in 1:length(counts_30)){
  average_30 = c(average_30, 
                 sum(counts_30[[i]]$count))
}

average_56 = NULL
for (i in 1:length(counts_56)){
  average_56 = c(average_56, 
                 sum(counts_56[[i]]$count))
}

average_82 = NULL
for (i in 1:length(counts_82)){
  average_82 = c(average_82, 
                 sum(counts_82[[i]]$count))
}

data = data.frame('mean_insertions' = c(average_16, 
                                        average_30, 
                                        average_56, 
                                        average_82), 
                  'age' = c(rep('16 weeks', 
                                 length(counts_16)), 
                             rep("30 weeks", 
                                 length(counts_30)), 
                             rep("56 weeks", 
                                 length(counts_56)), 
                             rep("82 weeks", 
                                 length(counts_82))))

ggplot(data) + 
  geom_violin(aes(x = age, y = mean_insertions)) + 
  geom_boxplot(aes(x = age, y = mean_insertions), width = 0.1) + 
  ylim(c(0, 3)) +
  theme_bw()

ggsave("subset_mean_insertions.png", 
       height = 4, width = 5)

### files -----
files_16 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_16_weeks/C03/kidney_multiome/effect_size", 
                      full.names = FALSE)
files_30 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_30_weeks/L12/kidney_multiome/effect_size", 
                      full.names = FALSE)
files_56 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_56_weeks/Q17/kidney_multiome/effect_size", 
                      full.names = FALSE)
files_82 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_82_weeks/U21/kidney_multiome/effect_size", 
                      full.names = FALSE)
files = Reduce(intersect, list(files_16,
                               files_30,
                               files_56, 
                               files_82
                               ))
files = gtools::mixedsort(files)

logp_u21_multiome_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/logp_mean_counts.txt"))
logp_q17_multiome_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_56_weeks/Q17/kidney_multiome/logp_mean_counts_freezed.txt"))
logp_l12_multiome_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_30_weeks/L12/kidney_multiome/logp_mean_counts_freezed.txt"))
logp_c03_multiome_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_16_weeks/C03/kidney_multiome/logp_mean_counts_freezed.txt"))

### 16 weeks -----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_16_weeks/C03/kidney_multiome/effect_size_10")
binding_sites_C03 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_c03_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_c03_multiome_mean$labels_mean_counts))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)
# binding_sites_C03 = data.table::rbindlist(binding_sites_C03)
saveRDS(binding_sites_C03, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/age_16_weeks/C03/kidney_multiome/thresholded_binding_sites_10.rds')

### 30 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_30_weeks/L12/kidney_multiome/effect_size_10")
binding_sites_L12 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_l12_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_l12_multiome_mean$labels_mean_counts))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)
# binding_sites_L12 = data.table::rbindlist(binding_sites_L12)
saveRDS(binding_sites_L12, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/age_30_weeks/L12/kidney_multiome/thresholded_binding_sites_subset_10.rds')

### 56 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_56_weeks/Q17/kidney_multiome/effect_size_all")
binding_sites_Q17 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_q17_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_q17_multiome_mean$labels_mean_counts))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)

# binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17)
saveRDS(binding_sites_Q17, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/age_56_weeks/Q17/kidney_multiome/thresholded_binding_sites_all_10.rds')

### 82 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/age_82_weeks/U21/kidney_multiome/effect_size_all")
binding_sites_U21 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$total_coverage = binding_sites$mean_coverage*binding_sites$width
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_u21_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_u21_multiome_mean$labels_mean_counts))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)

# binding_sites_U21 = data.table::rbindlist(binding_sites_U21)
saveRDS(binding_sites_U21, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/age_82_weeks/U21/kidney_multiome/thresholded_binding_sites_all_10.rds')

### visualizations ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/")

binding_sites_C03 = readRDS('age_16_weeks/C03/kidney_multiome/thresholded_binding_sites.rds')
binding_sites_L12 = readRDS('age_30_weeks/L12/kidney_multiome/thresholded_binding_sites_subset.rds')
binding_sites_Q17 = readRDS('age_56_weeks/Q17/kidney_multiome/thresholded_binding_sites_subset.rds')
binding_sites_U21 = readRDS('age_82_weeks/U21/kidney_multiome/thresholded_binding_sites_subset.rds')

widths = seq(15, 200, 10)
binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill = T)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill = T)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill = T)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill = T)

count_footprinting_C03 = NULL
count_footprinting_L12 = NULL
count_footprinting_Q17 = NULL
count_footprinting_U21 = NULL
widths = seq(15, 200, 10)
for (j in widths){
  count_footprinting_C03 = c(count_footprinting_C03, 
                             dim(binding_sites_C03[binding_sites_C03$width_effect_size < j + 8 & 
                                                     binding_sites_C03$width_effect_size >= j - 7, ])[1])
  count_footprinting_L12 = c(count_footprinting_L12, 
                             dim(binding_sites_L12[binding_sites_L12$width_effect_size < j + 8 & 
                                                     binding_sites_L12$width_effect_size >= j - 7])[1])
  count_footprinting_Q17 = c(count_footprinting_Q17, 
                             dim(binding_sites_Q17[binding_sites_Q17$width_effect_size < j + 8 & 
                                                     binding_sites_Q17$width_effect_size >= j - 7])[1])
  count_footprinting_U21 = c(count_footprinting_U21, 
                             dim(binding_sites_U21[binding_sites_U21$width_effect_size < j + 8 & 
                                                     binding_sites_U21$width_effect_size >= j - 7])[1])
  
}
peaks_width = regionsBed$end - regionsBed$start + 1

results = data.frame('width' = widths, 
                              'weeks_16' = count_footprinting_C03*1000/sum(peaks_width), 
                              'weeks_30' = count_footprinting_L12*1000/sum(peaks_width), 
                              'weeks_56' = count_footprinting_Q17*1000/sum(peaks_width), 
                              'weeks_82' = count_footprinting_U21*1000/sum(peaks_width)
                     )
results = pivot_longer(results, 
                       cols = c('weeks_16', 
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82'
                                ),  
                       names_to = 'age',         
                       values_to = 'count')
results$age = factor(results$age, 
                     levels = c('weeks_16',
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82')
                     )

p1 = ggplot(results) + 
  geom_point(aes(x = width, y = count, color = age), 
             size = 1.5) + 
  geom_line(aes(x = width, y = count, color = age), 
            linewidth = 1.5) + 
  ylab("number of footprintings per kb") + 
  scale_x_continuous(limits = c(10, 200), 
                     breaks = c(10, seq(50, 200, 50))) + 
  labs(colour = "Age") + 
  scale_color_manual(values = c('#D3E671',
                                '#A3D8FF', 
                                '#FAA301',
                                '#EB8C93'),
                     labels = c("weeks_16" = "16 weeks", 
                                "weeks_30" = "30 weeks", 
                                "weeks_56" = "56 weeks", 
                                "weeks_82" = "82 weeks")) +
  theme_Publication()
p1
ggsave("age_num_comparison_subset_95.png", p1, 
       height = 5, width = 5)
ggsave("age_num_comparison_subset_95.pdf", p1, 
       height = 5, width = 5)

counts_sum = data.frame('counts' = c(sum(count_footprinting_C03[1:6]), 
                                  sum(count_footprinting_L12[1:6]), 
                                  sum(count_footprinting_Q17[1:6]), 
                                  sum(count_footprinting_U21[1:6])), 
                        'age' = c('16 weeks', '30 weeks', 
                                  '56 weeks', '82 weeks'), 
                        'se' = c(sqrt(sum(count_footprinting_C03[1:6])), 
                                 sqrt(sum(count_footprinting_L12[1:6])), 
                                 sqrt(sum(count_footprinting_Q17[1:6])), 
                                 sqrt(sum(count_footprinting_U21[1:6]))))

ggplot(counts_sum, aes(x = age, 
                       y = log10(counts))) +
  geom_bar(aes(x = age, 
               y = log10(counts), 
               fill = age), 
           stat = "identity", alpha = 0.9) +
  labs(y = "log count of TF footprints") +
  geom_errorbar(aes(ymax = log10(counts + 1.96*se), 
                    ymin = log10(counts - 1.96*se)),
                position=position_dodge(width=0.9), 
                width = 0.2, linewidth = 0.5) + 
  scale_fill_manual(values = c('#D3E671',
                               '#A3D8FF', 
                               '#FAA301',
                               '#E58C94')) +
  theme_Publication()

ggsave('counts_comparison_footprint_<70bps.pdf', 
       height = 4, width = 6)

### binding sites ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus")

binding_sites_C03 = readRDS('age_16_weeks/C03/kidney_multiome/thresholded_binding_sites.rds')
binding_sites_L12 = readRDS('age_30_weeks/L12/kidney_multiome/thresholded_binding_sites_subset.rds')
binding_sites_Q17 = readRDS('age_56_weeks/Q17/kidney_multiome/thresholded_binding_sites_subset.rds')
binding_sites_U21 = readRDS('age_82_weeks/U21/kidney_multiome/thresholded_binding_sites_subset.rds')

binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill=TRUE)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill=TRUE)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill=TRUE)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill=TRUE)

binding_sites_C03$start = binding_sites_C03$position_effect_size - binding_sites_C03$width_effect_size/2
binding_sites_C03$end = binding_sites_C03$position_effect_size + binding_sites_C03$width_effect_size/2
binding_sites_L12$start = binding_sites_L12$position_effect_size - binding_sites_L12$width_effect_size/2
binding_sites_L12$end = binding_sites_L12$position_effect_size + binding_sites_L12$width_effect_size/2
binding_sites_Q17$start = binding_sites_Q17$position_effect_size - binding_sites_Q17$width_effect_size/2
binding_sites_Q17$end = binding_sites_Q17$position_effect_size + binding_sites_Q17$width_effect_size/2
binding_sites_U21$start = binding_sites_U21$position_effect_size - binding_sites_U21$width_effect_size/2
binding_sites_U21$end = binding_sites_U21$position_effect_size + binding_sites_U21$width_effect_size/2

### hypomethylation -----
rat_clock_by_human <- readRDS('../clock_hypomethylation_rat.rds')

regionsBed = read.table("peaks_region.bed")
peaks = makeGRangesFromDataFrame(regionsBed)
aging_loci = unique(findOverlapPairs(rat_clock_by_human, 
                                     peaks)@first)
binding_sites_C03_ranges = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_L12_ranges = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_Q17_ranges = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_U21_ranges = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)

overlaps_loci_c03 = findOverlaps(aging_loci, binding_sites_C03_ranges)
overlaps_loci_l12 = findOverlaps(aging_loci, binding_sites_L12_ranges)
overlaps_loci_q17 = findOverlaps(aging_loci, binding_sites_Q17_ranges)
overlaps_loci_u21 = findOverlaps(aging_loci, binding_sites_U21_ranges)

loci_c03 = unique(aging_loci[overlaps_loci_c03@from])
loci_l12 = unique(aging_loci[overlaps_loci_l12@from])
loci_q17 = unique(aging_loci[overlaps_loci_q17@from])
loci_u21 = unique(aging_loci[overlaps_loci_u21@from])
c(length(loci_c03), length(loci_l12), 
  length(loci_q17), length(loci_u21))

### nucleosome
binding_sites_C03_nucle = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size > 120, ])
binding_sites_L12_nucle = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size > 120, ])
binding_sites_Q17_nucle = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size > 120, ])
binding_sites_U21_nucle = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size > 120, ])

overlaps_loci_c03_nucle = findOverlaps(aging_loci, binding_sites_C03_nucle)
overlaps_loci_l12_nucle = findOverlaps(aging_loci, binding_sites_L12_nucle)
overlaps_loci_q17_nucle = findOverlaps(aging_loci, binding_sites_Q17_nucle)
overlaps_loci_u21_nucle = findOverlaps(aging_loci, binding_sites_U21_nucle)
c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
  length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from))

loci_c03_nucle = aging_loci[overlaps_loci_c03_nucle@from]
loci_l12_nucle = aging_loci[overlaps_loci_l12_nucle@from]
loci_q17_nucle = aging_loci[overlaps_loci_q17_nucle@from]
loci_u21_nucle = aging_loci[overlaps_loci_u21_nucle@from]

length(findOverlaps(loci_c03_nucle, loci_l12_nucle)@from)/length(loci_c03_nucle)
length(findOverlaps(loci_l12_nucle, loci_q17_nucle)@from)/length(loci_l12_nucle)
length(findOverlaps(loci_q17_nucle, loci_u21_nucle)@from)/length(loci_q17_nucle)

overlaps_loci = data.frame('age' = c("weeks_16", 
                                     "weeks_30", 
                                     "weeks_56", 
                                     "weeks_82"), 
                           'TF' = c(length(overlaps_loci_c03@from), length(overlaps_loci_l12@from), 
                                    length(overlaps_loci_q17@from), length(overlaps_loci_u21@from)), 
                           'Nucleosome' = c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
                                            length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from)))
overlaps_loci = overlaps_loci %>%
  pivot_longer(cols = c(TF, Nucleosome), 
               values_to = 'count', 
               names_to = 'Footprinting')
overlaps_loci$Footprinting = factor(overlaps_loci$Footprinting, 
                                    levels = c("TF", "Nucleosome"))
p2 = ggplot(overlaps_loci) + 
  geom_bar(aes(x = Footprinting, y = count, fill = age), 
           stat="identity", color = 'white',
           position="dodge", alpha = 0.8) + 
  ylab("counts of overlapped age-associated loci") + 
  scale_fill_manual(values =  c('#D3E671',
                                '#A3D8FF', 
                                '#FAA301',
                                '#E58C94'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  theme_Publication()
p2
ggsave("counts_overlapped_age-associated_hypomethylation_loci_95_effect_size_70.png", 
       p2, 
       height = 5, width = 5.5)
ggsave("counts_overlapped_age-associated_hypomethylation_loci_95_effect_size_70.pdf", 
       height = 5, width = 5.5)

### epitrace clockDML -----
mouse_clock_by_human <- readRDS('../mouse_clock_lifted_human_to_mouse_mm10.rds')
rat_clock_by_human <- easyLift::easyLiftOver(mouse_clock_by_human, 
                                             map = '../mm10ToRn7.over.chain')
rat_clock_by_human = data.frame(rat_clock_by_human)
reference = read.csv("~/nzhanglab/data/ParkerWilson_DEFND/refernece/GRCr8_GCF_036323735.1_sequence_info.csv")
reference = reference[c('refseqAccession', 'sequenceName')]
rat_clock_by_human = left_join(rat_clock_by_human, reference, 
                               by = join_by(seqnames == sequenceName))
rat_clock_by_human = rat_clock_by_human %>%
  dplyr::select(refseqAccession, start, end)
colnames(rat_clock_by_human)[1] = 'chr'
rat_clock_by_human = makeGRangesFromDataFrame(rat_clock_by_human)

regionsBed = read.table("peaks_region.bed")
peaks = makeGRangesFromDataFrame(regionsBed)
aging_loci = unique(findOverlapPairs(rat_clock_by_human, 
                                     peaks)@first)

binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill=TRUE)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill=TRUE)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill=TRUE)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill=TRUE)

binding_sites_C03$start = binding_sites_C03$position_effect_size - binding_sites_C03$width_effect_size/2
binding_sites_C03$end = binding_sites_C03$position_effect_size + binding_sites_C03$width_effect_size/2
binding_sites_L12$start = binding_sites_L12$position_effect_size - binding_sites_L12$width_effect_size/2
binding_sites_L12$end = binding_sites_L12$position_effect_size + binding_sites_L12$width_effect_size/2
binding_sites_Q17$start = binding_sites_Q17$position_effect_size - binding_sites_Q17$width_effect_size/2
binding_sites_Q17$end = binding_sites_Q17$position_effect_size + binding_sites_Q17$width_effect_size/2
binding_sites_U21$start = binding_sites_U21$position_effect_size - binding_sites_U21$width_effect_size/2
binding_sites_U21$end = binding_sites_U21$position_effect_size + binding_sites_U21$width_effect_size/2
binding_sites_C03_ranges = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_L12_ranges = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_Q17_ranges = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_U21_ranges = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)

overlaps_loci_c03 = findOverlaps(aging_loci, binding_sites_C03_ranges)
overlaps_loci_l12 = findOverlaps(aging_loci, binding_sites_L12_ranges)
overlaps_loci_q17 = findOverlaps(aging_loci, binding_sites_Q17_ranges)
overlaps_loci_u21 = findOverlaps(aging_loci, binding_sites_U21_ranges)

loci_c03 = unique(aging_loci[overlaps_loci_c03@from])
loci_l12 = unique(aging_loci[overlaps_loci_l12@from])
loci_q17 = unique(aging_loci[overlaps_loci_q17@from])
loci_u21 = unique(aging_loci[overlaps_loci_u21@from])
c(length(loci_c03), length(loci_l12), 
  length(loci_q17), length(loci_u21))

### nucleosome
binding_sites_C03_nucle = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size > 120, ])
binding_sites_L12_nucle = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size > 120, ])
binding_sites_Q17_nucle = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size > 120, ])
binding_sites_U21_nucle = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size > 120, ])

overlaps_loci_c03_nucle = findOverlaps(aging_loci, binding_sites_C03_nucle)
overlaps_loci_l12_nucle = findOverlaps(aging_loci, binding_sites_L12_nucle)
overlaps_loci_q17_nucle = findOverlaps(aging_loci, binding_sites_Q17_nucle)
overlaps_loci_u21_nucle = findOverlaps(aging_loci, binding_sites_U21_nucle)
c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
  length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from))

loci_c03_nucle = aging_loci[overlaps_loci_c03_nucle@from]
loci_l12_nucle = aging_loci[overlaps_loci_l12_nucle@from]
loci_q17_nucle = aging_loci[overlaps_loci_q17_nucle@from]
loci_u21_nucle = aging_loci[overlaps_loci_u21_nucle@from]

length(findOverlaps(loci_c03_nucle, loci_l12_nucle)@from)/length(loci_c03_nucle)
length(findOverlaps(loci_l12_nucle, loci_q17_nucle)@from)/length(loci_l12_nucle)
length(findOverlaps(loci_q17_nucle, loci_u21_nucle)@from)/length(loci_q17_nucle)

overlaps_loci = data.frame('age' = c("weeks_16", 
                                     "weeks_30", 
                                     "weeks_56", 
                                     "weeks_82"), 
                           'TF' = c(length(overlaps_loci_c03@from), length(overlaps_loci_l12@from), 
                                    length(overlaps_loci_q17@from), length(overlaps_loci_u21@from)), 
                           'Nucleosome' = c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
                                            length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from)))
overlaps_loci = overlaps_loci %>%
  pivot_longer(cols = c(TF, Nucleosome), 
               values_to = 'count', 
               names_to = 'Footprinting')
overlaps_loci$Footprinting = factor(overlaps_loci$Footprinting, 
                                    levels = c("TF", "Nucleosome"))
p2 = ggplot(overlaps_loci) + 
  geom_bar(aes(x = Footprinting, y = count, fill = age), 
           stat="identity", color = 'white',
           position="dodge", alpha = 0.8) + 
  ylab("counts of overlapped age-associated loci") + 
  scale_fill_manual(values =  c('#D3E671',
                                '#A3D8FF', 
                                '#FAA301',
                                '#E58C94'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  theme_Publication()
p2
ggsave("counts_overlapped_age-associated_loci_95_effect_size_70.png", 
       p2, 
       height = 5, width = 5.5)
ggsave("counts_overlapped_age-associated_loci_95_effect_size_70.pdf", 
       height = 5, width = 5.5)

### PRINT loci -----
## c03 ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/age_16_weeks/C03/kidney_multiome")
files = list.files("chunkedTFBSResults", full.names = TRUE)
chunkResults = NULL
for (i in 1:length(files)){
  chunkResults = c(chunkResults, 
                   readRDS(files[i]))
}
saveRDS(chunkResults, 
        'chunkedTFBSResults/chunk.rds')
chunkResults = readRDS("chunkedTFBSResults/chunk.rds")
binding_sites_multiome_c03_PRINT = pbmcapply::pbmclapply(seq_along(chunkResults), 
                                                         function(i){
                                                           regionTFBS = chunkResults[[i]]
                                                           sites = regionTFBS$sites
                                                           sites$TFBSScores = regionTFBS$TFBSScores
                                                           sites = sites[sites$TFBSScores >= 0.5]
                                                           reduced_sites = GenomicRanges::reduce(sites)
                                                           binding_sites = data.frame(reduced_sites)
                                                           binding_sites = binding_sites %>%
                                                             dplyr::select(seqnames, start, end, width)
                                                           return(binding_sites)
                                                         }, mc.cores = 2)

saveRDS(binding_sites_multiome_c03_PRINT, 
        "binding_sites_multiome_c03_PRINT.rds")
## l12 ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/age_30_weeks/L12/kidney_multiome")
files = list.files("chunkedTFBSResults", full.names = TRUE)
chunkResults = NULL
for (i in 1:length(files)){
  chunkResults = c(chunkResults, 
                   readRDS(files[i]))
}
saveRDS(chunkResults, 
        'chunkedTFBSResults/chunk.rds')
# chunkResults = readRDS("chunkedTFBSResults/chunk.rds")
binding_sites_multiome_l12_PRINT = pbmcapply::pbmclapply(seq_along(chunkResults), 
                                                         function(i){
                                                           regionTFBS = chunkResults[[i]]
                                                           sites = regionTFBS$sites
                                                           sites$TFBSScores = regionTFBS$TFBSScores
                                                           sites = sites[sites$TFBSScores >= 0.5]
                                                           reduced_sites = GenomicRanges::reduce(sites)
                                                           binding_sites = data.frame(reduced_sites)
                                                           binding_sites = binding_sites %>%
                                                             dplyr::select(seqnames, start, end, width)
                                                           return(binding_sites)
                                                         }, mc.cores = 4)

saveRDS(binding_sites_multiome_l12_PRINT, 
        "binding_sites_multiome_l12_PRINT.rds")
## q17 ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/age_30_weeks/L12/kidney_multiome")
files = list.files("chunkedTFBSResults", full.names = TRUE)
chunkResults = NULL
for (i in 1:length(files)){
  chunkResults = c(chunkResults, 
                   readRDS(files[i]))
}
saveRDS(chunkResults, 
        'chunkedTFBSResults/chunk.rds')
# chunkResults = readRDS("chunkedTFBSResults/chunk.rds")
binding_sites_multiome_q17_PRINT = pbmcapply::pbmclapply(seq_along(chunkResults), 
                                                         function(i){
                                                           regionTFBS = chunkResults[[i]]
                                                           sites = regionTFBS$sites
                                                           sites$TFBSScores = regionTFBS$TFBSScores
                                                           sites = sites[sites$TFBSScores >= 0.5]
                                                           reduced_sites = GenomicRanges::reduce(sites)
                                                           binding_sites = data.frame(reduced_sites)
                                                           binding_sites = binding_sites %>%
                                                             dplyr::select(seqnames, start, end, width)
                                                           return(binding_sites)
                                                         }, mc.cores = 4)

saveRDS(binding_sites_multiome_q17_PRINT, 
        "binding_sites_multiome_q17_PRINT.rds")
## u21 ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/age_82_weeks/U21/kidney_multiome")
files = list.files("chunkedTFBSResults", full.names = TRUE)
chunkResults = NULL
for (i in 1:length(files)){
  chunkResults = c(chunkResults, 
                   readRDS(files[i]))
}
saveRDS(chunkResults, 
        'chunkedTFBSResults/chunk.rds')
# chunkResults = readRDS("chunkedTFBSResults/chunk.rds")
binding_sites_multiome_u21_PRINT = pbmcapply::pbmclapply(seq_along(chunkResults), 
                                                         function(i){
                                                           regionTFBS = chunkResults[[i]]
                                                           sites = regionTFBS$sites
                                                           sites$TFBSScores = regionTFBS$TFBSScores
                                                           sites = sites[sites$TFBSScores >= 0.5]
                                                           reduced_sites = GenomicRanges::reduce(sites)
                                                           binding_sites = data.frame(reduced_sites)
                                                           binding_sites = binding_sites %>%
                                                             dplyr::select(seqnames, start, end, width)
                                                           return(binding_sites)
                                                          }, mc.cores = 4)


saveRDS(binding_sites_multiome_u21_PRINT, 
        "binding_sites_multiome_u21_PRINT.rds")

### PRINT visualization ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus")
binding_sites_multiome_c03_PRINT = readRDS('age_16_weeks/C03/kidney_multiome/binding_sites_multiome_c03_PRINT.rds')
binding_sites_multiome_l12_PRINT = readRDS('age_30_weeks/L12/kidney_multiome/binding_sites_multiome_l12_PRINT.rds')
binding_sites_multiome_q17_PRINT = readRDS('age_56_weeks/Q17/kidney_multiome/binding_sites_multiome_q17_PRINT.rds')
binding_sites_multiome_u21_PRINT = readRDS('age_82_weeks/U21/kidney_multiome/binding_sites_multiome_u21_PRINT.rds')

binding_sites_multiome_c03_PRINT = data.table::rbindlist(binding_sites_multiome_c03_PRINT)
binding_sites_multiome_l12_PRINT = data.table::rbindlist(binding_sites_multiome_l12_PRINT)
binding_sites_multiome_q17_PRINT = data.table::rbindlist(binding_sites_multiome_q17_PRINT)
binding_sites_multiome_u21_PRINT = data.table::rbindlist(binding_sites_multiome_u21_PRINT)

widths = seq(15, 200, 10)
count_footprinting_C03_PRINT = NULL
count_footprinting_L12_PRINT = NULL
count_footprinting_Q17_PRINT = NULL
count_footprinting_U21_PRINT = NULL

for (j in widths){
  count_footprinting_C03_PRINT = c(count_footprinting_C03_PRINT, 
                                   dim(binding_sites_multiome_c03_PRINT[binding_sites_multiome_c03_PRINT$width < j + 5 & 
                                                                          binding_sites_multiome_c03_PRINT$width >= j - 5, ])[1])
  count_footprinting_L12_PRINT = c(count_footprinting_L12_PRINT, 
                                   dim(binding_sites_multiome_l12_PRINT[binding_sites_multiome_l12_PRINT$width < j + 5 & 
                                                                          binding_sites_multiome_l12_PRINT$width >= j - 5])[1])
  count_footprinting_Q17_PRINT = c(count_footprinting_Q17_PRINT, 
                                   dim(binding_sites_multiome_q17_PRINT[binding_sites_multiome_q17_PRINT$width < j + 5 & 
                                                                          binding_sites_multiome_q17_PRINT$width >= j - 5])[1])
  count_footprinting_U21_PRINT = c(count_footprinting_U21_PRINT, 
                                   dim(binding_sites_multiome_u21_PRINT[binding_sites_multiome_u21_PRINT$width < j + 5 & 
                                                                          binding_sites_multiome_u21_PRINT$width >= j - 5])[1])
  
}
peaks_width = regionsBed$end - regionsBed$start + 1

results = data.frame('width' = widths, 
                     'weeks_16' = count_footprinting_C03_PRINT*1000/sum(peaks_width), 
                     'weeks_30' = count_footprinting_L12_PRINT*1000/sum(peaks_width), 
                     'weeks_56' = count_footprinting_Q17_PRINT*1000/sum(peaks_width), 
                     'weeks_82' = count_footprinting_U21_PRINT*1000/sum(peaks_width)
)
results = pivot_longer(results, 
                       cols = c('weeks_16', 
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82'
                       ),  
                       names_to = 'age',         
                       values_to = 'count')
results$age = factor(results$age, 
                     levels = c('weeks_16',
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82')
)

p1 = ggplot(results) + 
  geom_point(aes(x = width, y = count, color = age), 
             size = 1.5) + 
  geom_line(aes(x = width, y = count, color = age), 
            linewidth = 1.5) + 
  ylab("footprintings number per kb") + 
  scale_x_continuous(limits = c(10, 200), 
                     breaks = c(10, seq(50, 200, 50))) + 
  labs(colour = "Age") + 
  scale_color_manual(values = c('#D3E671',
                                '#A3D8FF', 
                                '#FAA301',
                                '#F37199'), 
                     labels = c("weeks_16" = "16 weeks", 
                                "weeks_30" = "30 weeks", 
                                "weeks_56" = "56 weeks", 
                                "weeks_82" = "82 weeks")) +
  theme_Publication()
p1
ggsave("age_num_comparison_subset_PRINT.pdf", 
       height = 5, width = 5)
ggsave("age_num_comparison_subset_PRINT.png", 
       height = 5, width = 5)

mouse_clock_by_human <- readRDS('../mouse_clock_lifted_human_to_mouse_mm10.rds')
rat_clock_by_human <- easyLift::easyLiftOver(mouse_clock_by_human, 
                                             map = '../mm10ToRn7.over.chain')
rat_clock_by_human = data.frame(rat_clock_by_human)
reference = read.csv("~/nzhanglab/data/ParkerWilson_DEFND/refernece/GRCr8_GCF_036323735.1_sequence_info.csv")
reference = reference[c('refseqAccession', 'sequenceName')]
rat_clock_by_human = left_join(rat_clock_by_human, reference, 
                               by = join_by(seqnames == sequenceName))
rat_clock_by_human = rat_clock_by_human %>%
  dplyr::select(refseqAccession, start, end)
colnames(rat_clock_by_human)[1] = 'chr'
rat_clock_by_human = makeGRangesFromDataFrame(rat_clock_by_human)

regionsBed = read.table("peaks_region.bed")
peaks = makeGRangesFromDataFrame(regionsBed)
aging_loci = unique(findOverlapPairs(rat_clock_by_human, 
                                     peaks)@first)

binding_sites_multiome_c03_PRINT_ranges = makeGRangesFromDataFrame(binding_sites_multiome_c03_PRINT[binding_sites_multiome_c03_PRINT$width <= 80, ])
binding_sites_multiome_l12_PRINT_ranges = makeGRangesFromDataFrame(binding_sites_multiome_l12_PRINT[binding_sites_multiome_l12_PRINT$width <= 80, ])
binding_sites_multiome_q17_PRINT_ranges = makeGRangesFromDataFrame(binding_sites_multiome_q17_PRINT[binding_sites_multiome_q17_PRINT$width <= 80, ])
binding_sites_multiome_u21_PRINT_ranges = makeGRangesFromDataFrame(binding_sites_multiome_u21_PRINT[binding_sites_multiome_u21_PRINT$width <= 80, ])

overlaps_loci_c03_PRINT = findOverlaps(aging_loci, binding_sites_multiome_c03_PRINT_ranges)
overlaps_loci_l12_PRINT = findOverlaps(aging_loci, binding_sites_multiome_l12_PRINT_ranges)
overlaps_loci_q17_PRINT = findOverlaps(aging_loci, binding_sites_multiome_q17_PRINT_ranges)
overlaps_loci_u21_PRINT = findOverlaps(aging_loci, binding_sites_multiome_u21_PRINT_ranges)

loci_c03_PRINT = aging_loci[overlaps_loci_c03_PRINT@from]
loci_l12_PRINT = aging_loci[overlaps_loci_l12_PRINT@from]
loci_q17_PRINT = aging_loci[overlaps_loci_q17_PRINT@from]
loci_u21_PRINT = aging_loci[overlaps_loci_u21_PRINT@from]

length(loci_c03_PRINT)
length(loci_l12_PRINT)
length(loci_q17_PRINT)
length(loci_u21_PRINT)

binding_sites_multiome_c03_PRINT_nucle_ranges = makeGRangesFromDataFrame(binding_sites_multiome_c03_PRINT[binding_sites_multiome_c03_PRINT$width >= 120 & 
                                                                                                            binding_sites_multiome_c03_PRINT$width <= 200, ])
binding_sites_multiome_l12_PRINT_nucle_ranges = makeGRangesFromDataFrame(binding_sites_multiome_l12_PRINT[binding_sites_multiome_l12_PRINT$width >= 120 & 
                                                                                                            binding_sites_multiome_l12_PRINT$width <= 200, ])
binding_sites_multiome_q17_PRINT_nucle_ranges = makeGRangesFromDataFrame(binding_sites_multiome_q17_PRINT[binding_sites_multiome_q17_PRINT$width >= 120 & 
                                                                                                            binding_sites_multiome_q17_PRINT$width <= 200, ])
binding_sites_multiome_u21_PRINT_nucle_ranges = makeGRangesFromDataFrame(binding_sites_multiome_u21_PRINT[binding_sites_multiome_u21_PRINT$width >= 120 & 
                                                                                                            binding_sites_multiome_u21_PRINT$width <= 200, ])

overlaps_loci_c03_PRINT_nucle = findOverlaps(aging_loci, binding_sites_multiome_c03_PRINT_nucle_ranges)
overlaps_loci_l12_PRINT_nucle = findOverlaps(aging_loci, binding_sites_multiome_l12_PRINT_nucle_ranges)
overlaps_loci_q17_PRINT_nucle = findOverlaps(aging_loci, binding_sites_multiome_q17_PRINT_nucle_ranges)
overlaps_loci_u21_PRINT_nucle = findOverlaps(aging_loci, binding_sites_multiome_u21_PRINT_nucle_ranges)

loci_c03_PRINT_nucle = aging_loci[overlaps_loci_c03_PRINT_nucle@from]
loci_l12_PRINT_nucle = aging_loci[overlaps_loci_l12_PRINT_nucle@from]
loci_q17_PRINT_nucle = aging_loci[overlaps_loci_q17_PRINT_nucle@from]
loci_u21_PRINT_nucle = aging_loci[overlaps_loci_u21_PRINT_nucle@from]

length(loci_c03_PRINT_nucle)
length(loci_l12_PRINT_nucle)
length(loci_q17_PRINT_nucle)
length(loci_u21_PRINT_nucle)

overlaps_loci_PRINT = data.frame('age' = c("weeks_16", 
                                           "weeks_30", 
                                           "weeks_56", 
                                           "weeks_82"), 
                                 'TF' = c(length(overlaps_loci_c03_PRINT@from), length(overlaps_loci_l12_PRINT@from), 
                                          length(overlaps_loci_q17_PRINT@from), length(overlaps_loci_u21_PRINT@from)), 
                                 'Nucleosome' = c(length(loci_c03_PRINT_nucle), length(loci_l12_PRINT_nucle), 
                                                  length(loci_q17_PRINT_nucle), length(loci_u21_PRINT_nucle)))
overlaps_loci_PRINT = overlaps_loci_PRINT %>%
  pivot_longer(cols = c(TF, Nucleosome), 
               values_to = 'count', 
               names_to = 'Footprinting')
# greenBlue = c('#e0f4cc','#ccebc5','#a8ddc1',"#a1cbc1")
overlaps_loci_PRINT$Footprinting = factor(overlaps_loci_PRINT$Footprinting, 
                                          levels = c("TF", "Nucleosome"))
p1 = ggplot(overlaps_loci_PRINT) + 
  geom_bar(aes(x = Footprinting, y = count, fill = age), 
           stat="identity", 
           position="dodge", alpha = 0.9) + 
  # scale_fill_manual(values = greenBlue) + 
  ylab("counts of overlapped age-associated loci by PRINT") + 
  scale_fill_manual(values = c('#D3E671',
                               '#A3D8FF', 
                               '#FAA301',
                               '#F37199'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  theme_Publication()
p1
ggsave("counts_overlapped_age-associated_loci_PRINT.pdf", 
       height = 5, width = 6)
ggsave("counts_overlapped_age-associated_loci.png", 
       height = 5, width = 6)

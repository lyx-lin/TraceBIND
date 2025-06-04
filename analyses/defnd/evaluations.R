library(tidyverse)
library(ggpubr)
library(GenomicRanges)
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

### age 82 weeks -----
setwd(paste0("~/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/"))

source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("peaks_region.txt")
colnames(regionsBed) = c('chr', 'start', 'end')
### NBinomial ----
files_multiome_u21 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome_subset_effect_size_freezed2", 
                                full.names = FALSE)
files_defnd_u21 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/U21/kidney_defnd_subset_effect_size_freezed2", 
                             full.names = FALSE)
files_defnd_v22 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/V22/kidney_defnd_effect_size_freezed", 
                             full.names = FALSE)
files = Reduce(intersect, list(files_multiome_u21,
                               files_defnd_u21,
                               files_defnd_v22))
files = gtools::mixedsort(files)
index = as.numeric(sapply(files, 
                          function(x) strsplit(x, '_')[[1]][1]))
logp_u21_multiome_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/logp_mean_counts_freezed.txt"))
logp_u21_defnd_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_defnd/logp_mean_counts_freezed.txt"))
logp_v22_defnd_mean = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/V22/kidney_defnd/logp_mean_counts_freezed.txt"))

setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome_subset_effect_size_freezed")
binding_sites_multiome_u21 = pbmcapply::pbmclapply(1:length(files), 
                                function(i) {
                                  binding_sites = read.table(files[i], header = TRUE) 
                                  j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                  if (dim(binding_sites)[1] > 0){
                                    binding_sites$threshold = apply(binding_sites, 1, function(x)
                                      logp_u21_multiome_mean$threshold[which.min(abs(x['coverage'] - 
                                                                                       as.numeric(as.character(logp_u21_multiome_mean$labels_mean_counts))))])
                                    binding_sites = binding_sites[binding_sites$p_value >= 
                                                                    pmax(-log10(0.05), binding_sites$threshold), ]                                  
                                  }
                                  binding_sites$chr = rep(regionsBed$chr[j], dim(binding_sites)[1])
                                  
                                  binding_sites
                                }, 
                                mc.cores = 2
                                )
binding_sites_multiome_u21 = data.table::rbindlist(binding_sites_multiome_u21, 
                                                   fill=TRUE)
saveRDS(binding_sites_multiome_u21, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/thresholded_binding_sites_subset.rds')

setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/U21/kidney_defnd_subset_effect_size_freezed")
binding_sites_defnd_u21 = pbmcapply::pbmclapply(1:length(files), 
                                                   function(i) {
                                                     binding_sites = read.table(files[i], header = TRUE) 
                                                     j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                                     if (dim(binding_sites)[1] > 0){
                                                       binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                         logp_u21_defnd_mean$smoothed_threshold[which.min(abs(x['coverage'] - 
                                                                                                                as.numeric(as.character(logp_u21_defnd_mean$labels_mean_counts))))])
                                                       binding_sites = binding_sites[binding_sites$p_value >= 
                                                                                       pmax(-log10(0.05), binding_sites$threshold), ]                                                     
                                                     }
                                                     binding_sites$chr = rep(regionsBed$chr[j], 
                                                                             dim(binding_sites)[1])
                                                     binding_sites
                                                   }, 
                                                
                                                   mc.cores = 2
)
binding_sites_defnd_u21 = data.table::rbindlist(binding_sites_defnd_u21, 
                                                fill = T)
saveRDS(binding_sites_defnd_u21, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_defnd/thresholded_binding_sites_subset.rds')

setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/age_82_weeks/V22/kidney_defnd_effect_size_freezed")
binding_sites_defnd_v22 = pbmcapply::pbmclapply(1:length(files), 
                                                function(i) {
                                                  binding_sites = read.table(files[i], header = TRUE) 
                                                  j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                                  if (dim(binding_sites)[1] > 0){
                                                    binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                      logp_v22_defnd_mean$smoothed_threshold[which.min(abs(x['coverage'] - 
                                                                                                             as.numeric(as.character(logp_v22_defnd_mean$labels_mean_counts))))])
                                                    binding_sites = binding_sites[binding_sites$p_value >= 
                                                                                    pmax(-log10(0.05), binding_sites$threshold), ]
                                                  }
                                                  binding_sites$chr = rep(regionsBed$chr[j], dim(binding_sites)[1])
                                                  
                                                  binding_sites
                                                }, 
                                                mc.cores = 2
)
binding_sites_defnd_v22 = data.table::rbindlist(binding_sites_defnd_v22, fill = T)

saveRDS(binding_sites_defnd_v22, 
        '~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/V22/kidney_defnd/thresholded_binding_sites.rds')

## counts ----
binding_sites_defnd_v22 = readRDS('~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/V22/kidney_defnd/thresholded_binding_sites.rds')
binding_sites_defnd_u21 = readRDS('~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_defnd/thresholded_binding_sites_subset.rds')
binding_sites_multiome_u21 = readRDS('~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/thresholded_binding_sites_subset.rds')

setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks")
widths = seq(15, 200, 10)

count_footprinting_multiome_u21 = NULL
count_footprinting_defnd_u21 = NULL
count_footprinting_defnd_v22 = NULL

for (j in widths){
  count_footprinting_multiome_u21 = c(count_footprinting_multiome_u21, 
                                      dim(binding_sites_multiome_u21[binding_sites_multiome_u21$width_effect_size < j + 8 & 
                                                                       binding_sites_multiome_u21$width_effect_size >= j - 3, ])[1])
  count_footprinting_defnd_u21 = c(count_footprinting_defnd_u21, 
                                   dim(binding_sites_defnd_u21[binding_sites_defnd_u21$width_effect_size < j + 8 & 
                                                                 binding_sites_defnd_u21$width_effect_size >= j - 3, ])[1])
  count_footprinting_defnd_v22 = c(count_footprinting_defnd_v22, 
                                   dim(binding_sites_defnd_v22[binding_sites_defnd_v22$width_effect_size < j + 8 & 
                                                                 binding_sites_defnd_v22$width_effect_size >= j - 3, ])[1])
}

colnames(regionsBed) = c('chr', 'start', 'end')
peaks_width = sum(regionsBed$end - regionsBed$start + 1)
index = as.numeric(sapply(files, function(x) strsplit(x, "_")[[1]][1]))

results = data.frame('width' = widths, 
                     'No_depletion' = count_footprinting_multiome_u21*1000/sum(regionsBed$end[index] - regionsBed$start[index] + 1), 
                     'Weak_depletion' = count_footprinting_defnd_u21*1000/sum(regionsBed$end[index] - regionsBed$start[index] + 1), 
                     'Strong_depletion' = count_footprinting_defnd_v22*1000/sum(regionsBed$end[index] - regionsBed$start[index] + 1))
results = pivot_longer(results, 
                       cols = c("No_depletion", 
                                'Weak_depletion', 
                                'Strong_depletion'),  
                       names_to = 'type',         
                       values_to = 'count')
results$type = factor(results$type, 
                      levels = c("Strong_depletion", 
                                 "Weak_depletion", 
                                 "No_depletion"))
p1 = ggplot(results) + 
  geom_point(aes(x = width, y = count, color = type), 
             size = 1, alpha = 1) + 
  geom_line(aes(x = width, y = count, color = type), 
             linewidth = 1, alpha = 1) + 
  ylab("number of TraceBind footprintings per kb") + 
  scale_x_continuous(limits = c(10, 200), 
                     breaks = c(10, seq(50, 200, 50))) + 
  labs(colour = "Depletion Level") + 
  scale_color_manual(values = c("#859F3D", 
                                "#C6D332", 
                                "#FFB0B0"), 
                     labels = c("No_depletion" = "No depletion", 
                                "Weak_depletion" = "Weak depletion", 
                                "Strong_depletion" = "Strong depletion")) +
  theme_Publication()

p1
ggsave("NBinomial_counts_comparison_subset.png", 
       p1, height = 6, width = 6)
ggsave("NBinomial_counts_comparison_subset.pdf", 
       p1, height = 6, width = 6)
### enrichment of nucleosome  ----
binding_sites_multiome_u21_nucleosome = binding_sites_multiome_u21[binding_sites_multiome_u21$width >= 120, ]
binding_sites_multiome_u21_nucleosome$start = binding_sites_multiome_u21_nucleosome$position - 
  binding_sites_multiome_u21_nucleosome$width/2
binding_sites_multiome_u21_nucleosome$end = binding_sites_multiome_u21_nucleosome$position + 
  binding_sites_multiome_u21_nucleosome$width/2

binding_sites_multiome_u21_nucleosome_granges = makeGRangesFromDataFrame(binding_sites_multiome_u21_nucleosome, 
                                                                         keep.extra.columns=TRUE)

binding_sites_defnd_u21_nucleosome = binding_sites_defnd_u21[binding_sites_defnd_u21$width >= 120, ]
binding_sites_defnd_u21_nucleosome$start = binding_sites_defnd_u21_nucleosome$position - 
  binding_sites_defnd_u21_nucleosome$width/2
binding_sites_defnd_u21_nucleosome$end = binding_sites_defnd_u21_nucleosome$position + 
  binding_sites_defnd_u21_nucleosome$width/2

binding_sites_defnd_u21_nucleosome_granges = makeGRangesFromDataFrame(binding_sites_defnd_u21_nucleosome, 
                                                                      keep.extra.columns=TRUE)

binding_sites_defnd_v22_nucleosome = binding_sites_defnd_v22[binding_sites_defnd_v22$width >= 120, ]
binding_sites_defnd_v22_nucleosome$start = binding_sites_defnd_v22_nucleosome$position - 
  binding_sites_defnd_v22_nucleosome$width/2
binding_sites_defnd_v22_nucleosome$end = binding_sites_defnd_v22_nucleosome$position + 
  binding_sites_defnd_v22_nucleosome$width/2

binding_sites_defnd_v22_nucleosome_granges = makeGRangesFromDataFrame(binding_sites_defnd_v22_nucleosome, 
                                                                      keep.extra.columns=TRUE)

overlaps_multiome_u21 = findOverlapPairs(binding_sites_multiome_u21_nucleosome_granges, 
                                         binding_sites_defnd_u21_nucleosome_granges)@first
multiome_u21_nucle = overlaps_multiome_u21[-findOverlaps(overlaps_multiome_u21, binding_sites_defnd_v22_nucleosome_granges)@from]
defnd_u21_nucle = binding_sites_defnd_u21_nucleosome_granges[findOverlaps(multiome_u21_nucle, binding_sites_defnd_u21_nucleosome_granges)@to]
multiome_u21_nucle = multiome_u21_nucle[findOverlaps(multiome_u21_nucle, binding_sites_defnd_u21_nucleosome_granges)@from]
range = multiome_u21_nucle
peaks = makeGRangesFromDataFrame(regionsBed)

overlaps_defnd_u21 = findOverlapPairs(binding_sites_multiome_u21_nucleosome_granges, 
                                      binding_sites_defnd_u21_nucleosome_granges)@second
overlaps = findOverlaps(overlaps_multiome_u21, binding_sites_defnd_v22_nucleosome_granges)
overlaps_multiome_u21_nucleosome = overlaps_multiome_u21[queryHits(overlaps)]
overlaps_defnd_u21_nucleosome = overlaps_defnd_u21[queryHits(overlaps)]
overlaps_defnd_v22_nucleosome = binding_sites_defnd_v22_nucleosome_granges[subjectHits(overlaps)]

data = data.frame('enrichment' = c(overlaps_defnd_v22_nucleosome$effect_size, 
                                   overlaps_defnd_u21_nucleosome$effect_size, 
                                   overlaps_multiome_u21_nucleosome$effect_size), 
                  'depletion' = c(rep("Strong_depletion", length(overlaps_defnd_v22_nucleosome)), 
                                  rep("Weak_depletion", length(overlaps_defnd_u21_nucleosome)), 
                                  rep("No_depletion", length(overlaps_multiome_u21_nucleosome))))

data$depletion = factor(data$depletion, 
                        levels = c("Strong_depletion", 
                                   "Weak_depletion", 
                                   "No_depletion"))
logfc_data <- data %>%
  group_by(depletion) %>%
  summarise(mean_enrichment = mean(enrichment)) %>%
  arrange(depletion) 
strong_weak = -log2(logfc_data$mean_enrichment[logfc_data$depletion == 'Strong_depletion']/logfc_data$mean_enrichment[logfc_data$depletion == 'Weak_depletion'])
strong_no = -log2(logfc_data$mean_enrichment[logfc_data$depletion == 'Strong_depletion']/logfc_data$mean_enrichment[logfc_data$depletion == 'No_depletion'])
weak_no = -log2(logfc_data$mean_enrichment[logfc_data$depletion == 'Weak_depletion']/logfc_data$mean_enrichment[logfc_data$depletion == 'No_depletion'])

p3 = ggplot(data, aes(x = depletion, y = enrichment)) +
  # geom_violin() +
  geom_boxplot(aes(fill = depletion), show.legend = FALSE) +
  ylim(c(0, 2.5)) + 
  labs(x = "Depletion Level",
       y = "TraceBind Effect size") +
  scale_x_discrete(labels = c("No_depletion" = "No depletion", 
                              "Weak_depletion" = "Weak depletion", 
                              "Strong_depletion" = "Strong depletion")) +
  scale_fill_manual(values = c("#859F3D", 
                               "#C6D332", 
                               "#FFB0B0"), 
                     labels = c("No_depletion" = "No depletion", 
                                "Weak_depletion" = "Weak depletion", 
                                "Strong_depletion" = "Strong depletion")) +
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.signif", 
                     comparisons = list(c("Weak_depletion", "Strong_depletion"),
                                        c("No_depletion", "Weak_depletion"),
                                        c("No_depletion", "Strong_depletion")),
                     method.args = list(alternative = "greater"),
                     tip.length = 0.01, 
                     label.y = c(1.2, 1.4, 1.6)
                     ) + 
  annotate('text', x = 1.6, y = 2.5, size = 4, 
           label = "One-sided Paired Wilcoxon Test") + 
  theme_Publication()

p3
ggsave("NBinomial_enrichment_comparison_120_subset.png", 
       p3, height = 3, width = 4)
ggsave("NBinomial_enrichment_comparison_120_subset.pdf", 
       p3, height = 3, width = 4)
ks.test(overlaps_defnd_u21_nucleosome$effect_size, 
        overlaps_multiome_u21_nucleosome$effect_size, 
        alternative = 'greater')

ks.test(overlaps_defnd_v22_nucleosome$effect_size, 
        overlaps_multiome_u21_nucleosome$effect_size, 
        alternative = 'greater')

ks.test(overlaps_defnd_v22_nucleosome$effect_size, 
        overlaps_defnd_u21_nucleosome$effect_size, 
        alternative = 'greater')


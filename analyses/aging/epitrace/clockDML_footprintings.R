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

setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all"))
source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')

## subpeak seurat =====
library(EpiTrace)
library(Seurat)
library(Signac)
library(tidyverse)
library(GenomicFeatures)

setwd("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")

seuratobject_fp = readRDS("seuratobject_merged_split_peaks_no_smoothed_PT_thresholded_95_effect_size.rds")
seuratobject_fp = seuratobject_fp %>%
  subset(celltype %in% c("PT", "PT-MT"))

binding_sites_all = readRDS("thresholded_binding_sites_PT.rds")
mouse_clock_by_human <- readRDS('../mouse_clock_lifted_human_to_mouse_mm10.rds')
rat_clock_by_human <- easyLift::easyLiftOver(mouse_clock_by_human, map = '../mm10ToRn7.over.chain')
rat_clock_by_human = data.frame(rat_clock_by_human)
reference = read.csv("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/refernece/GRCr8_GCF_036323735.1_sequence_info.csv")
reference = reference[c('refseqAccession', 'sequenceName')]
rat_clock_by_human = left_join(rat_clock_by_human, reference, 
                               by = join_by(seqnames == sequenceName))
rat_clock_by_human = rat_clock_by_human %>%
  dplyr::select(refseqAccession, start, end)
colnames(rat_clock_by_human)[1] = 'chr'
rat_clock_by_human = makeGRangesFromDataFrame(rat_clock_by_human)

binding_sites_all = data.table::rbindlist(binding_sites_all, fill = T)
binding_sites_all = binding_sites_all[binding_sites_all$p_value >= -log10(0.05), ]
binding_sites_all$start = binding_sites_all$position_effect_size - binding_sites_all$width_effect_size/2
binding_sites_all$end = binding_sites_all$position_effect_size + binding_sites_all$width_effect_size/2
binding_sites_all_ranges = makeGRangesFromDataFrame(binding_sites_all[binding_sites_all$width_effect_size <= 120 & 
                                                                        binding_sites_all$width_effect_size >= 8,  ], 
                                                    keep.extra.columns = T)
overlaps_loci = findOverlaps(rat_clock_by_human, binding_sites_all_ranges)
locis = rat_clock_by_human[overlaps_loci@from]

locis = unique(locis)

locis = data.frame(locis)
locis = left_join(locis, reference, 
                  by = join_by(seqnames == refseqAccession))
locis = locis %>%
  dplyr::select(sequenceName, start, end)
colnames(locis)[1] = 'chr'
locis = makeGRangesFromDataFrame(locis)

# aging_loci = data.frame(aging_loci)
# aging_loci = left_join(aging_loci, reference, 
#                        by = join_by(seqnames == refseqAccession))
# aging_loci = aging_loci %>%
#   dplyr::select(sequenceName, start, end)
# colnames(aging_loci)[1] = 'chr'
# aging_loci = makeGRangesFromDataFrame(aging_loci)

split_peaks = data.frame(granges(seuratobject_fp))
split_peaks = left_join(split_peaks, reference, 
                        by = join_by(seqnames == refseqAccession))
split_peaks = split_peaks %>%
  dplyr::select(sequenceName, start, end)
colnames(split_peaks)[1] = 'chr'

init_gr <- Init_Peakset(split_peaks)
matrix = seuratobject_fp[['ATAC']]@counts
rownames(matrix) = paste(split_peaks$chr, split_peaks$start, 
                         split_peaks$end, sep = '-')
init_mm <- Init_Matrix(peakname = paste(split_peaks$chr, split_peaks$start, 
                                        split_peaks$end, sep = '-'),
                       cellname = colnames(seuratobject_fp),
                       matrix = matrix)

### subpeaks & footprinting -----
epitrace_obj_age_conv_human_clock_fp <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                matrix = init_mm,
                                                                ref_genome = 'rn7',
                                                                clock_gr = locis,
                                                                iterative_time = 1,
                                                                min.cutoff = 0,
                                                                non_standard_clock = T,
                                                                ncore_lim = 1,
                                                                mean_error_limit = 0.1)
epitrace_obj_age_conv_human_clock_fp@meta.data <- merge(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                                                        seuratobject_fp@meta.data, 
                                                        by=0)

logfc_data_fp_init <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_16']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30'])
weeks30_56_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56'])
weeks56_82_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_82'])
c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)

logfc_data_fp_iter <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_16']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30'])
weeks30_56_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56'])
weeks56_82_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_82'])
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

epitrace_obj_age_conv_human_clock_fp@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_fp@meta.data$dataset)

epitrace_fp_initial = ggplot(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                             aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.9) + 
  geom_boxplot(width=0.15) + 
  # geom_violin(data = metadata_initials_random, 
  #             aes(x = dataset, 
  #                 y = EpiTraceAge_Clock_initial, 
  #                 fill = dataset), 
  #             alpha = 0.4, color = 'gray') + 
  scale_fill_manual(values =  c('#B7D433',
                                '#8ECEF1', 
                                '#FCAF17',
                                '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Epitrace Initial Age") + 
  xlab('Age') + 
  theme_Publication()

epitrace_fp_iteration = ggplot(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                               aes(x=dataset, y=EpiTraceAge_iterative)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.9) + 
  geom_boxplot(width=0.15) + 
  # geom_boxplot(data = epitrace_obj_age_conv_human_clock_random@meta.data, 
  #              aes(x = dataset, 
  #                  y = EpiTraceAge_iterative, 
  #                  fill = dataset), alpha = 0.2, 
  #              width=0.15) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Epitrace Iterative Age") + 
  xlab('Age') + 
  theme_Publication()
epitrace_fp_initial
epitrace_fp_iteration

c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

ggsave("epitrace_fp_split_width_120_effect_size.pdf", 
       epitrace_fp_initial + epitrace_fp_iteration + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"), 
       height = 5, width = 9)
ggsave("epitrace_fp_initial_split_width_120_effect_size.pdf", 
       epitrace_fp_initial, 
       height = 5, width = 6)
saveRDS(epitrace_obj_age_conv_human_clock_fp, 
        "epitrace_obj_age_conv_human_clock_fp_split_no_smoothed_thresholded_95_effect_size.rds")

### subpeaks & all clockDML ----
epitrace_obj_age_conv_human_clock_all <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                 matrix = init_mm,
                                                                 ref_genome = 'rn7',
                                                                 clock_gr = aging_loci,
                                                                 iterative_time = 1,
                                                                 min.cutoff = 0,
                                                                 non_standard_clock = F,
                                                                 qualnum = 10,
                                                                 ncore_lim = 2,
                                                                 mean_error_limit = 0.1)

epitrace_obj_age_conv_human_clock_all@meta.data <- merge(epitrace_obj_age_conv_human_clock_all@meta.data, 
                                                         seuratobject_split@meta.data, 
                                                         by=0)
logfc_data_all_init <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_16']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30'])
weeks30_56_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56'])
weeks56_82_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_82'])
c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init)

logfc_data_all_iter <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_16']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30'])
weeks30_56_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56'])
weeks56_82_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_82'])
c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter)

epitrace_obj_age_conv_human_clock_all@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_all@meta.data$dataset)

epitrace_all_initial = ggplot(epitrace_obj_age_conv_human_clock_all@meta.data, 
                              aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
  geom_violin(scale='width',aes(fill=dataset)) + 
  geom_boxplot(width=0.2,
               fill='black',outlier.alpha = 0) + 
  ggpubr::stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 3, 
                     comparisons = list(c("weeks_16", "weeks_30"),
                                        c("weeks_30", "weeks_56"),
                                        c("weeks_56", "weeks_82")),
                     method.args = list(alternative = "less"),
                     tip.length = 0.03) + 
  theme_classic()

epitrace_all_iteration = ggplot(epitrace_obj_age_conv_human_clock_all@meta.data, 
                                aes(x=dataset, y=EpiTraceAge_iterative)) + 
  geom_violin(scale='width',aes(fill=dataset)) + 
  geom_boxplot(width=0.2,fill='black',outlier.alpha = 0) + 
  ggpubr::stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 3, 
                     comparisons = list(c("weeks_16", "weeks_30"),
                                        c("weeks_30", "weeks_56"),
                                        c("weeks_56", "weeks_82")),
                     method.args = list(alternative = "less"),
                     tip.length = 0.03) + 
  theme_classic()
c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init)
c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter)

ggsave("epitrace_all_iteration_split.pdf", 
       epitrace_all_iteration, 
       height = 5, width = 8)
ggsave("epitrace_all_initial_split.pdf", 
       epitrace_all_initial, 
       height = 5, width = 8)
saveRDS(epitrace_obj_age_conv_human_clock_all, 
        "epitrace_obj_age_conv_human_clock_all_split.rds")

### peak seurat =====
library(EpiTrace)
library(Seurat)
library(Signac)

setwd("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")

seuratobject_whole = readRDS("seuratobject_whole_peaks.rds")
seuratobject_whole = seuratobject_whole %>%
  subset(celltype %in% c("PT", "PT-MT"))

locis = unique(locis)
locis = data.frame(locis)
locis = left_join(locis, reference, 
                  by = join_by(seqnames == refseqAccession))
locis = locis %>%
  dplyr::select(sequenceName, start, end)
colnames(locis)[1] = 'chr'
locis = makeGRangesFromDataFrame(locis)

aging_loci = data.frame(aging_loci)
aging_loci = left_join(aging_loci, reference, 
                       by = join_by(seqnames == refseqAccession))
aging_loci = aging_loci %>%
  dplyr::select(sequenceName, start, end)
colnames(aging_loci)[1] = 'chr'
aging_loci = makeGRangesFromDataFrame(aging_loci)

peaks = data.frame(granges(seuratobject_whole))
peaks = left_join(peaks, reference, 
                        by = join_by(seqnames == refseqAccession))
peaks = peaks %>%
  dplyr::select(sequenceName, start, end)
colnames(peaks)[1] = 'chr'

init_gr <- Init_Peakset(peaks)
matrix = seuratobject_whole[['ATAC']]@counts
rownames(matrix) = paste(peaks$chr, peaks$start, 
                         peaks$end, sep = '-')
init_mm <- Init_Matrix(peakname = paste(peaks$chr, peaks$start, 
                                        peaks$end, sep = '-'),
                       cellname = colnames(seuratobject_whole),
                       matrix = matrix)
### all peaks & footprinting ----
epitrace_obj_age_conv_human_clock_fp <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                matrix = init_mm,
                                                                ref_genome = 'rn7',
                                                                clock_gr = locis,
                                                                iterative_time = 1,
                                                                min.cutoff = 0,
                                                                non_standard_clock = T,
                                                                # qualnum = 10,
                                                                ncore_lim = 1,
                                                                mean_error_limit = 0.1)
epitrace_obj_age_conv_human_clock_fp@meta.data <- merge(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                                                        seuratobject_whole@meta.data, 
                                                        by=0)

logfc_data_fp_init <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_16']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30'])
weeks30_56_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56'])
weeks56_82_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_82'])
c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)

logfc_data_fp_iter <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_16']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30'])
weeks30_56_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56'])
weeks56_82_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_82'])
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

epitrace_obj_age_conv_human_clock_fp@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_fp@meta.data$dataset)

epitrace_fp_initial = ggplot(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                             aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
  geom_violin(scale='width',aes(fill=dataset), alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) +   
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 3, 
                     comparisons = list(c("weeks_16", "weeks_30"),
                                        c("weeks_30", "weeks_56"),
                                        c("weeks_56", "weeks_82")),
                     method.args = list(alternative = "less"),
                     tip.length = 0.03) + 
  theme_classic()

epitrace_fp_iteration = ggplot(epitrace_obj_age_conv_human_clock_fp@meta.data, 
                               aes(x=dataset, y=EpiTraceAge_iterative)) + 
  geom_violin(scale='width',aes(fill=dataset), alpha = 0.8) + 
  geom_boxplot(width=0.2, outlier.alpha = 0) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 3, 
                     comparisons = list(c("weeks_16", "weeks_30"),
                                        c("weeks_30", "weeks_56"),
                                        c("weeks_56", "weeks_82")),
                     method.args = list(alternative = "less"),
                     tip.length = 0.03) + 
  theme_Publication()
# epitrace_fp_initial
# epitrace_fp_iteration
c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

ggsave("epitrace_whole_fp.pdf", 
       epitrace_fp_initial + epitrace_fp_iteration + 
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"), 
       height = 5, width = 10)

ggsave("epitrace_fp_iteration_whole_peak_width_120.pdf", 
       epitrace_fp_iteration, 
       height = 5, width = 5)
ggsave("epitrace_fp_initial_whole_peak_width_120.pdf", 
       epitrace_fp_initial, 
       height = 5, width = 5)
# saveRDS(epitrace_obj_age_conv_human_clock_fp, 
#         "epitrace_obj_age_conv_human_clock_fp.rds")

### all peaks & all clockdmls ----
epitrace_obj_age_conv_human_clock_all <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                 matrix = init_mm,
                                                                 ref_genome = 'rn7',
                                                                 clock_gr = aging_loci,
                                                                 iterative_time = 1,
                                                                 min.cutoff = 0,
                                                                 non_standard_clock = F,
                                                                 qualnum = 10,
                                                                 ncore_lim = 1,
                                                                 mean_error_limit = 0.1)

epitrace_obj_age_conv_human_clock_all@meta.data <- merge(epitrace_obj_age_conv_human_clock_all@meta.data, 
                                                         seuratobject_split@meta.data, 
                                                         by=0)
logfc_data_all_init <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_16']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30'])
weeks30_56_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56'])
weeks56_82_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_82'])
c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init)

logfc_data_all_iter <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_16']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30'])
weeks30_56_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56'])
weeks56_82_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_82'])
c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter)

epitrace_obj_age_conv_human_clock_all@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_all@meta.data$dataset)

epitrace_all_initial = ggplot(epitrace_obj_age_conv_human_clock_all@meta.data, 
                              aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                     labels = c("weeks_16" = "16 weeks", 
                                "weeks_30" = "30 weeks", 
                                "weeks_56" = "56 weeks", 
                                "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks"
                             )) + 
  ylab("Epitrace Initial Age") + 
  xlab('Age') + 
  theme_Publication()

epitrace_all_iteration = ggplot(epitrace_obj_age_conv_human_clock_all@meta.data, 
                                aes(x=dataset, y=EpiTraceAge_iterative)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Epitrace Iterative Age") + 
  xlab('Age') + 
  theme_Publication()
c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init)
c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter)

epitrace_all_initial + epitrace_all_iteration + 
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("epitrace_all.pdf", 
       epitrace_all_initial + epitrace_all_iteration + 
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"), 
       height = 5, width = 10)
ggsave("epitrace_all_iteration.pdf", 
       epitrace_all_iteration, 
       height = 3, width = 4)
### all peaks & random clockdmls -----
index = findOverlaps(aging_loci, locis)@from
setwd("epitrace_downsample_all")
for (k in 1:30){
  set.seed(k)
  index = findOverlaps(aging_loci, locis)@from
  epitrace_obj_age_conv_human_clock_random <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                      matrix = init_mm,
                                                                      ref_genome = 'rn7',
                                                                      clock_gr = aging_loci[sample(setdiff(1:length(aging_loci), index), 530)],
                                                                      iterative_time = 1,
                                                                      min.cutoff = 0,
                                                                      non_standard_clock = T,
                                                                      qualnum = 10,
                                                                      ncore_lim = 2,
                                                                      mean_error_limit = 0.1)
  epitrace_obj_age_conv_human_clock_random@meta.data <- merge(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                                              seuratobject@meta.data, 
                                                              by=0)
  epitrace_obj_age_conv_human_clock_random@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_random@meta.data$dataset)
  
  logfc_data_random_initial <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_Clock_initial, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_16']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30'])
  weeks30_56_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56'])
  weeks56_82_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_82'])
  c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial)
  
  epitrace_random_initial = ggplot(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                   aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
    geom_violin(scale='width',aes(fill=dataset)) + 
    annotate('text', x = 1.5, y = 1.1, 
             label = round(weeks16_30_random_initial, 2)) + 
    annotate('text', x = 2.5, y = 1.1, 
             label = round(weeks30_56_random_initial, 2)) + 
    annotate('text', x = 3.5, y = 1.1, 
             label = round(weeks56_82_random_initial, 2)) + 
    geom_boxplot(width=0.2, fill='white',
                 outlier.alpha = 0) + 
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               label = "p.signif", 
                               size = 3, 
                               comparisons = list(c("weeks_16", "weeks_30"),
                                                  c("weeks_30", "weeks_56"),
                                                  c("weeks_56", "weeks_82")),
                               method.args = list(alternative = "less"),
                               tip.length = 0.03) + 
    theme_classic()
  ggsave(paste0("epitrace_random_all_initial_", k, ".pdf"), 
         epitrace_random_initial, 
         height = 5, width = 8)
  logfc_data_random_iterative <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_iterative, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_16']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30'])
  weeks30_56_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56'])
  weeks56_82_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_82'])
  c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative)
  
  epitrace_random_iteration = ggplot(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                     aes(x=dataset, y=EpiTraceAge_iterative)) + 
    geom_violin(scale='width',aes(fill=dataset)) + 
    geom_boxplot(width=0.2, fill='black',outlier.alpha = 0) + 
    annotate('text', x = 1.5, y = 1.1, 
             label = round(weeks16_30_random_iterative, 2)) + 
    annotate('text', x = 2.5, y = 1.1, 
             label = round(weeks30_56_random_iterative, 2)) + 
    annotate('text', x = 3.5, y = 1.1, 
             label = round(weeks56_82_random_iterative, 2)) + 
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               label = "p.signif", 
                               size = 3, 
                               comparisons = list(c("weeks_16", "weeks_30"),
                                                  c("weeks_30", "weeks_56"),
                                                  c("weeks_56", "weeks_82")),
                               method.args = list(alternative = "less"),
                               tip.length = 0.03) + 
    theme_classic()
  ggsave(paste0("epitrace_random_all_iteration_", k, ".pdf"), 
         epitrace_random_iteration, 
         height = 5, width = 8)
}

### subpeaks & random clockdmls -----
index = findOverlaps(aging_loci, locis)@from
setwd("epitrace_downsample_splitup")
split_peaks = data.frame(granges(seuratobject_fp))
split_peaks = left_join(split_peaks, reference, 
                        by = join_by(seqnames == refseqAccession))
split_peaks = split_peaks %>%
  dplyr::select(sequenceName, start, end)
colnames(split_peaks)[1] = 'chr'

init_gr <- Init_Peakset(split_peaks)
matrix = seuratobject_fp[['ATAC']]@counts
rownames(matrix) = paste(split_peaks$chr, split_peaks$start, 
                         split_peaks$end, sep = '-')
init_mm <- Init_Matrix(peakname = paste(split_peaks$chr, split_peaks$start, 
                                        split_peaks$end, sep = '-'),
                       cellname = colnames(seuratobject_fp),
                       matrix = matrix)

for (k in 30:30){
  set.seed(k)
  index = findOverlaps(aging_loci, locis)@from
  epitrace_obj_age_conv_human_clock_random <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                      matrix = init_mm,
                                                                      ref_genome = 'rn7',
                                                                      clock_gr = aging_loci[sample(setdiff(1:length(aging_loci), index), 530)],
                                                                      iterative_time = 1,
                                                                      min.cutoff = 0,
                                                                      non_standard_clock = T,
                                                                      qualnum = 10,
                                                                      ncore_lim = 2,
                                                                      mean_error_limit = 0.1)
  epitrace_obj_age_conv_human_clock_random@meta.data <- merge(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                                              seuratobject_fp@meta.data, 
                                                              by=0)
  epitrace_obj_age_conv_human_clock_random@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_random@meta.data$dataset)
  
  logfc_data_random_initial <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_Clock_initial, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_16']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30'])
  weeks30_56_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56'])
  weeks56_82_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_82'])
  c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial)
  
  epitrace_random_initial = ggplot(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                   aes(x=dataset, y=EpiTraceAge_Clock_initial)) + 
    geom_violin(scale='width',aes(fill=dataset)) + 
    annotate('text', x = 1.5, y = 1.1, 
             label = round(weeks16_30_random_initial, 2)) + 
    annotate('text', x = 2.5, y = 1.1, 
             label = round(weeks30_56_random_initial, 2)) + 
    annotate('text', x = 3.5, y = 1.1, 
             label = round(weeks56_82_random_initial, 2)) + 
    geom_boxplot(width=0.2, fill='white',
                 outlier.alpha = 0) + 
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               label = "p.signif", 
                               size = 3, 
                               comparisons = list(c("weeks_16", "weeks_30"),
                                                  c("weeks_30", "weeks_56"),
                                                  c("weeks_56", "weeks_82")),
                               method.args = list(alternative = "less"),
                               tip.length = 0.03) + 
    theme_classic()
  ggsave(paste0("epitrace_random_splitup_initial_", k, ".pdf"), 
         epitrace_random_initial, 
         height = 5, width = 8)
  logfc_data_random_iterative <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_iterative, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_16']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30'])
  weeks30_56_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56'])
  weeks56_82_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_82'])
  c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative)
  
  epitrace_random_iteration = ggplot(epitrace_obj_age_conv_human_clock_random@meta.data, 
                                     aes(x=dataset, y=EpiTraceAge_iterative)) + 
    geom_violin(scale='width',aes(fill=dataset)) + 
    geom_boxplot(width=0.2, fill='black',outlier.alpha = 0) + 
    annotate('text', x = 1.5, y = 1.1, 
             label = round(weeks16_30_random_iterative, 2)) + 
    annotate('text', x = 2.5, y = 1.1, 
             label = round(weeks30_56_random_iterative, 2)) + 
    annotate('text', x = 3.5, y = 1.1, 
             label = round(weeks56_82_random_iterative, 2)) + 
    ggpubr::stat_compare_means(method = "wilcox.test", 
                               label = "p.signif", 
                               size = 3, 
                               comparisons = list(c("weeks_16", "weeks_30"),
                                                  c("weeks_30", "weeks_56"),
                                                  c("weeks_56", "weeks_82")),
                               method.args = list(alternative = "less"),
                               tip.length = 0.03) + 
    theme_classic()
  ggsave(paste0("epitrace_random_splitup_iteration_", k, ".pdf"), 
         epitrace_random_iteration, 
         height = 5, width = 8)
  saveRDS(epitrace_obj_age_conv_human_clock_random,
          paste0('epitrace_obj_age_conv_human_clock_random_splitup_', k, '.rds'))
}
## random whole v.s. footprinting comparison ----
epitrace_obj_age_conv_human_clock_fp = readRDS("epitrace_obj_age_conv_human_clock_fp_split_no_smoothed_thresholded.rds")
logfc_data_fp_init <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_16']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30'])
weeks30_56_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56'])
weeks56_82_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_82'])
c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)

logfc_data_fp_iter <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_16']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30'])
weeks30_56_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56'])
weeks56_82_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_82'])
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

epitrace_obj_age_conv_human_clock_fp@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_fp@meta.data$dataset)

initials_logfc_random = NULL
iterative_logfc_random = NULL
for (k in 1:15){
  epitrace_obj_age_conv_human_clock_random = readRDS(paste0("epitrace_downsample_all/epitrace_obj_age_conv_human_clock_random_all_", k, ".rds"))
  logfc_data_random_initial <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_Clock_initial, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_16']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30'])
  weeks30_56_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56'])
  weeks56_82_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_82'])
  initials_logfc_random = rbind(initials_logfc_random, 
                                c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial))
  
  logfc_data_random_iterative <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_iterative, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_16']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30'])
  weeks30_56_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56'])
  weeks56_82_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_82'])
  iterative_logfc_random = rbind(iterative_logfc_random, 
                                 c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative))
  
}
apply(initials_logfc_random, 2, mean)
# 0.13655961 0.09347811 0.07586573
apply(iterative_logfc_random, 2, mean)
# 0.1725314 0.1995416 0.1077074

metadata_initial_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_Clock_initial')]
metadata_iterative_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_iterative')]
for (k in 1:30){
  epitrace_obj_age_conv_human_clock_random = readRDS(paste0("epitrace_downsample_all/epitrace_obj_age_conv_human_clock_random_all_", k, ".rds"))
  initial = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_Clock_initial')]
  colnames(initial)[2] = k
  metadata_initial_random = merge(metadata_initial_random, 
                                  initial, 
                                  by='Row.names')
  iterative = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_iterative')]
  colnames(iterative)[2] = k
  metadata_iterative_random = merge(metadata_iterative_random, 
                                    iterative, 
                                    by='Row.names')
}

metadata_initial_random$random = apply(metadata_initial_random[, 4:33], 1, mean)
metadata_iterative_random$random = apply(metadata_iterative_random[, 4:33], 1, mean)

logfc_data_random_initial_mean <- metadata_initial_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_16']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30'])
weeks30_56_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56'])
weeks56_82_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_82'])

logfc_data_random_iterative_mean <- metadata_iterative_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_16']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30'])
weeks30_56_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56'])
weeks56_82_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_82'])

c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial)
# 0.13677981 0.09241749 0.06183347
c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative)
# 0.1583694 0.1903182 0.1008868
logfc_data_random_initial_mean$type = 'random'
epitrace_random_initial = ggplot(metadata_initial_random, 
                             aes(x=dataset, 
                                 y=random)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) +
  # geom_violin(data = metadata_initial_random, 
  #             aes(x = dataset, 
  #                 y = random), width = 0.6, 
  #             alpha = 0.4, color = 'gray') + 
  # geom_boxplot(data = metadata_initial_random, 
  #              aes(x = dataset, 
  #                  y = random), 
  #              width=0.15, alpha = 0.4, color = 'gray') +
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Predicted Age") + 
  xlab('Age') + 
  theme_Publication()
epitrace_random_initial
epitrace_random_iteration = ggplot(metadata_iteration_random, 
                               aes(x=dataset, y=random)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) + 
  # geom_violin(data = metadata_iterative_random, 
  #             aes(x = dataset, 
  #                 y = EpiTraceAge_iterative, 
  #                 fill = dataset), 
  #             alpha = 0.4, color = 'gray') + 
  # geom_boxplot(data = epitrace_obj_age_conv_human_clock_random@meta.data, 
  #              aes(x = dataset, 
  #                  y = EpiTraceAge_iterative, 
  #                  fill = dataset), alpha = 0.2, 
  #              width=0.15) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Epitrace Iterative Age") + 
  xlab('Age') + 
  theme_Publication()
epitrace_random_initial
epitrace_random_iteration

ggsave("epitrace_random_whole_width_120.pdf", 
       epitrace_random_initial + epitrace_random_iteration + 
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"), 
       height = 5, width = 10)

## random splitup v.s. footprinting comparison ----
epitrace_obj_age_conv_human_clock_fp = readRDS("epitrace_obj_age_conv_human_clock_fp_split_no_smoothed_thresholded.rds")
epitrace_obj_age_conv_human_clock_fp@meta.data$dataset = factor(epitrace_obj_age_conv_human_clock_fp@meta.data$dataset)

initials_logfc_random = NULL
iterative_logfc_random = NULL
for (k in 1:30){
  epitrace_obj_age_conv_human_clock_random = readRDS(paste0("epitrace_downsample_splitup/epitrace_obj_age_conv_human_clock_random_splitup_", k, ".rds"))
  logfc_data_random_initial <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_Clock_initial, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_16']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30'])
  weeks30_56_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_30']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56'])
  weeks56_82_random_initial = -log2(logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_56']/logfc_data_random_initial$mean_age[logfc_data_random_initial$dataset == 'weeks_82'])
  initials_logfc_random = rbind(initials_logfc_random, 
                                c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial))
  
  logfc_data_random_iterative <- epitrace_obj_age_conv_human_clock_random@meta.data %>%
    group_by(dataset) %>%
    summarise(mean_age = mean(EpiTraceAge_iterative, na.rm = T)) %>%
    arrange(dataset) 
  weeks16_30_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_16']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30'])
  weeks30_56_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_30']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56'])
  weeks56_82_random_iterative = -log2(logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_56']/logfc_data_random_iterative$mean_age[logfc_data_random_iterative$dataset == 'weeks_82'])
  iterative_logfc_random = rbind(iterative_logfc_random, 
                                 c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative))
  
}
apply(initials_logfc_random, 2, mean)
# 0.13655961 0.09347811 0.07586573
apply(iterative_logfc_random, 2, mean)
# 0.1725314 0.1995416 0.1077074

metadata_initial_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_Clock_initial')]
metadata_iterative_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_iterative')]
for (k in 1:30){
  epitrace_obj_age_conv_human_clock_random = readRDS(paste0("epitrace_downsample_splitup/epitrace_obj_age_conv_human_clock_random_splitup_", k, ".rds"))
  initial = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_Clock_initial')]
  colnames(initial)[2] = k
  metadata_initial_random = merge(metadata_initial_random, 
                                  initial, 
                                  by='Row.names')
  iterative = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_iterative')]
  colnames(iterative)[2] = k
  metadata_iterative_random = merge(metadata_iterative_random, 
                                    iterative, 
                                    by='Row.names')
}

metadata_initial_random$random = apply(metadata_initial_random[, 4:33], 1, mean)
metadata_iterative_random$random = apply(metadata_iterative_random[, 4:33], 1, mean)

logfc_data_random_initial_mean <- metadata_initial_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_16']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30'])
weeks30_56_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56'])
weeks56_82_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_82'])

logfc_data_random_iterative_mean <- metadata_iterative_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_16']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30'])
weeks30_56_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56'])
weeks56_82_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_82'])

c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial)
#  0.08202364 0.06860174 0.05898958
c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative)
# 0.1687433 0.4408847 0.1812196
logfc_data_random_initial_mean$type = 'random'
epitrace_random_initial = ggplot(metadata_initial_random, 
                                 aes(x=dataset, 
                                     y=random)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) +
  # geom_violin(data = metadata_initial_random, 
  #             aes(x = dataset, 
  #                 y = random), width = 0.6, 
  #             alpha = 0.4, color = 'gray') + 
  # geom_boxplot(data = metadata_initial_random, 
  #              aes(x = dataset, 
  #                  y = random), 
  #              width=0.15, alpha = 0.4, color = 'gray') +
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Predicted Age") + 
  xlab('Age') + 
  theme_Publication()
epitrace_random_initial
epitrace_random_iteration = ggplot(metadata_iterative_random, 
                                   aes(x=dataset, y=random)) + 
  geom_violin(scale='width',aes(fill=dataset), 
              alpha = 0.8) + 
  geom_boxplot(width=0.2, alpha = 0.8) + 
  # geom_violin(data = metadata_iterative_random, 
  #             aes(x = dataset, 
  #                 y = EpiTraceAge_iterative, 
  #                 fill = dataset), 
  #             alpha = 0.4, color = 'gray') + 
  # geom_boxplot(data = epitrace_obj_age_conv_human_clock_random@meta.data, 
  #              aes(x = dataset, 
  #                  y = EpiTraceAge_iterative, 
  #                  fill = dataset), alpha = 0.2, 
  #              width=0.15) + 
  scale_fill_manual(values = c('#B7D433',
                               '#8ECEF1', 
                               '#FCAF17',
                               '#FB8691'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             label = "p.signif", 
                             size = 3, 
                             comparisons = list(c("weeks_16", "weeks_30"),
                                                c("weeks_30", "weeks_56"),
                                                c("weeks_56", "weeks_82")),
                             method.args = list(alternative = "less"),
                             tip.length = 0.03) + 
  scale_x_discrete(labels= c("16 weeks", 
                             "30 weeks", 
                             "56 weeks", 
                             "82 weeks")) + 
  ylab("Epitrace Iterative Age") + 
  xlab('Age') + 
  theme_Publication()
epitrace_random_initial
epitrace_random_iteration

ggsave("epitrace_random_splitup_width_120.pdf", 
       epitrace_random_initial + epitrace_random_iteration + 
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom"), 
       height = 5, width = 10)

## random v.s. footprint v.s. whole ----
epitrace_obj_age_conv_human_clock_fp = readRDS("epitrace_obj_age_conv_human_clock_fp_split_no_smoothed_thresholded_95_effect_size.rds")
logfc_data_fp_init <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_16']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30'])
weeks30_56_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_30']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56'])
weeks56_82_fp_init = -log2(logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_56']/logfc_data_fp_init$mean_age[logfc_data_fp_init$dataset == 'weeks_82'])
c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init)

logfc_data_fp_iter <- epitrace_obj_age_conv_human_clock_fp@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_16']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30'])
weeks30_56_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_30']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56'])
weeks56_82_fp_iter = -log2(logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_56']/logfc_data_fp_iter$mean_age[logfc_data_fp_iter$dataset == 'weeks_82'])
c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter)

epitrace_obj_age_conv_human_clock_all = readRDS("epitrace_obj_age_conv_human_clock_all.rds")
logfc_data_all_init <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_Clock_initial)) %>%
  arrange(dataset) 
weeks16_30_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_16']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30'])
weeks30_56_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_30']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56'])
weeks56_82_all_init = -log2(logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_56']/logfc_data_all_init$mean_age[logfc_data_all_init$dataset == 'weeks_82'])
c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init)

logfc_data_all_iter <- epitrace_obj_age_conv_human_clock_all@meta.data %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(EpiTraceAge_iterative)) %>%
  arrange(dataset) 
weeks16_30_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_16']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30'])
weeks30_56_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_30']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56'])
weeks56_82_all_iter = -log2(logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_56']/logfc_data_all_iter$mean_age[logfc_data_all_iter$dataset == 'weeks_82'])
c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter)

metadata_initial_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_Clock_initial')]
metadata_iterative_random = epitrace_obj_age_conv_human_clock_fp@meta.data[c('Row.names', 'dataset', 'EpiTraceAge_iterative')]
for (k in 1:20){
  epitrace_obj_age_conv_human_clock_random = readRDS(paste0("epitrace_downsample_all/epitrace_obj_age_conv_human_clock_random_all_", k, ".rds"))
  initial = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_Clock_initial')]
  colnames(initial)[2] = k
  metadata_initial_random = merge(metadata_initial_random, 
                                  initial, 
                                  by='Row.names')
  iterative = epitrace_obj_age_conv_human_clock_random@meta.data[c('Row.names', 'EpiTraceAge_iterative')]
  colnames(iterative)[2] = k
  metadata_iterative_random = merge(metadata_iterative_random, 
                                    iterative, 
                                    by='Row.names')
}

metadata_initial_random$random = apply(metadata_initial_random[, 4:23], 1, mean)
metadata_iterative_random$random = apply(metadata_iterative_random[, 4:23], 1, mean)

logfc_data_random_initial_mean <- metadata_initial_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_16']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30'])
weeks30_56_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_30']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56'])
weeks56_82_random_initial = -log2(logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_56']/logfc_data_random_initial_mean$mean_age[logfc_data_random_initial_mean$dataset == 'weeks_82'])

logfc_data_random_iterative_mean <- metadata_iterative_random %>%
  group_by(dataset) %>%
  summarise(mean_age = mean(random, na.rm = T)) %>%
  arrange(dataset) 

weeks16_30_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_16']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30'])
weeks30_56_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_30']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56'])
weeks56_82_random_iterative = -log2(logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_56']/logfc_data_random_iterative_mean$mean_age[logfc_data_random_iterative_mean$dataset == 'weeks_82'])

c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial)
c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative)

logfc_inital = data.frame("random" = c(weeks16_30_random_initial, weeks30_56_random_initial, weeks56_82_random_initial), 
                          "all" = c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init), 
                          "footprint.split.up.peaks" = c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init), 
                          'round' = rep('initial', 3), 
                          'age' = c('16 v.s. 30', '30 v.s. 56', '56 v.s. 82')) 
logfc_iter = data.frame("random" = c(weeks16_30_random_iterative, weeks30_56_random_iterative, weeks56_82_random_iterative), 
                        "all" = c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter), 
                        "footprint.split.up.peaks" = c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter), 
                        'round' = rep('iterative', 3), 
                        'age' = c('16 v.s. 30', '30 v.s. 56', '56 v.s. 82')) 
logfc = rbind(logfc_inital, logfc_iter)
logfc = logfc %>%
  pivot_longer(cols = c(random, all, footprint.split.up.peaks), 
               names_to = "method", 
               values_to = 'logfc')

ggplot(logfc) + 
  geom_bar(aes(x = age, y = logfc, fill = method), 
           stat="identity", color = 'red',
           position="dodge", alpha = 0.8) + 
  facet_wrap(~round) + 
  scale_fill_manual(values = 
                        c("black", 
                          '#8d99a6',
                          'gray')) + 
  theme_Publication() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.margin=margin(t = 0, unit='cm'))

ggsave("logfc_ages.pdf", 
       height = 1.5, width = 4.5)
wilcox.test(c(weeks16_30_all_init, weeks30_56_all_init, weeks56_82_all_init), 
            c(weeks16_30_fp_init, weeks30_56_fp_init, weeks56_82_fp_init), 
            paired = T, alternative = 'greater')
wilcox.test(c(weeks16_30_all_iter, weeks30_56_all_iter, weeks56_82_all_iter), 
            c(weeks16_30_fp_iter, weeks30_56_fp_iter, weeks56_82_fp_iter), 
            paired = T, alternative = 'less')


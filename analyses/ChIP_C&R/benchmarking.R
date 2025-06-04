library(tidyverse)
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
            # legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            # legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

setwd("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking")
ctcf = readRDS("CTCF/data.rds")
nr3c1 = readRDS("NR3C1/human_kidney/Control6/data.rds")
nfat5_control = readRDS("NFAT5/Control/Control_6/data.rds")
nfat5_ckd = readRDS("NFAT5/CKD/CKD_4/data.rds")
hivep2_control = readRDS("HIVEP2/Control/Control_6/data.rds")
hivep2_ckd = readRDS("HIVEP2/CKD/data.rds")
global_ref = readRDS("../global_refer/Control_3/benchmarking_data.rds")
ctcf$transcription_factor = rep("CTCF", 5)
nr3c1$transcription_factor = rep("NR3C1", 5)
nfat5_control$transcription_factor = rep("NFAT5_Healthy", 5)
nfat5_ckd$transcription_factor = rep("NFAT5_CKD", 5)
hivep2_control$transcription_factor = rep("HIVEP2_Healthy", 5)
hivep2_ckd$transcription_factor = rep("HIVEP2_CKD", 5)
global_ref = global_ref[c(8, 7, 5, 4), c(1:4)]
global_ref$transcription_factor = rep("DNase_footprint", )
data = rbind(ctcf, nr3c1, 
             nfat5_control, nfat5_ckd, 
             hivep2_control, hivep2_ckd, global_ref)
data = data[data$method != 'Footprint-to-TF PRINT', ]
data$recovery = round(as.numeric(data$recovery), 3)
data$random_prob = round(as.numeric(data$random_prob), 2)

data$transcription_factor <- factor(data$transcription_factor, 
                                    levels = c('CTCF', "NR3C1", 
                                               "NFAT5_Healthy", "NFAT5_CKD", 
                                               "HIVEP2_Healthy", "HIVEP2_CKD", 
                                               "DNase_footprint"))

p1 = ggplot(data = data, aes(x = transcription_factor, y = recovery, 
                             fill = method)) + 
  geom_bar(aes(x = transcription_factor, y = recovery, 
               fill = method), 
           stat="identity", 
           position="dodge", alpha = 1) + 
  geom_text(aes(label=recovery, y = recovery), 
            position=position_dodge(width=0.9), vjust=0, size=5.3) +
  # geom_text(aes(label=random_prob, y = random_prob), 
  #           position=position_dodge(width=0.9), vjust=-0.5, size=5.3) +
  # geom_errorbar(aes(ymax = random_prob, ymin = random_prob, 
  #                   linetype = 'random'),
  #               position=position_dodge(width=0.9), width = 1) + 
  ylab("recovery rate") + 
  xlab("Transcription Factor") + 
  scale_fill_manual(values = brewer.pal(9, "Pastel1")[c(1, 2, 5, 4)]) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +  
  scale_linetype_manual("", values = c('random' = 'dashed')) + 
  theme(legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.text.x = element_text(color = "grey20", size = 16),
        axis.text.y = element_text(color = "grey20", size = 16),  
        axis.title.x = element_text(color = "grey20", size = 18),
        axis.title.y = element_text(color = "grey20", size = 18)) + 
  theme_Publication()
p1
ggsave("benchmark_result1.pdf", p1, 
       height = 5, width = 18)

p2 = ggplot(data = data, aes(x = transcription_factor, y = random_prob, 
                             fill = method)) + 
  geom_bar(aes(x = transcription_factor, y = random_prob, 
               fill = method), 
           stat="identity", 
           position="dodge", alpha = 1) + 
  geom_text(aes(label=random_prob, y = random_prob), 
            position=position_dodge(width=0.9), vjust=0, size=5.3) +
  # geom_text(aes(label=random_prob, y = random_prob), 
  #           position=position_dodge(width=0.9), vjust=-0.5, size=5.3) +
  # geom_errorbar(aes(ymax = random_prob, ymin = random_prob, 
  #                   linetype = 'random'),
  #               position=position_dodge(width=0.9), width = 1) + 
  ylab("random prob") + 
  xlab("Transcription Factor") + 
  scale_fill_manual(values = brewer.pal(9, "Pastel1")[c(1, 2, 5, 4)]) +
  scale_y_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.1)) +  
  # scale_linetype_manual("", values = c('random' = 'dashed')) + 
  theme(legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.text.x = element_text(color = "grey20", size = 16),
        axis.text.y = element_text(color = "grey20", size = 16),  
        axis.title.x = element_text(color = "grey20", size = 18),
        axis.title.y = element_text(color = "grey20", size = 18)) + 
  theme_Publication()
p2
ggsave("random_prob_result.pdf", p2, 
       height = 5, width = 18)
### DNase ----
global_ref = readRDS("../global_refer/Control_3/benchmarking_data.rds")
global_ref = global_ref[c(8, 7, 5, 4), ]
global_ref$total_num = as.numeric(global_ref$total_num)
global_ref$method = factor(global_ref$method, 
                           levels = c('NBinomial', 
                                      'seq2PRINT', 
                                      'HINT-ATAC', 
                                      'TOBIAS'))
ggplot(global_ref, 
       aes(x = method, y = total_num)) + 
  geom_bar(aes(fill = method), 
           stat="identity", 
           position="dodge", alpha = 1) + 
  geom_text(aes(label=total_num, y = total_num), 
            position=position_dodge(width=0.9), 
            vjust=0, size=4) +
  ylab("counts of identified footprints") + 
  scale_fill_manual(values = brewer.pal(9, "Pastel1")[c(1, 2, 5, 4)]) +
  theme(legend.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        axis.text.x = element_text(color = "grey20", size = 16),
        axis.text.y = element_text(color = "grey20", size = 16),  
        axis.title.x = element_text(color = "grey20", size = 18),
        axis.title.y = element_text(color = "grey20", size = 18)) + 
  theme_Publication()

ggsave('counts_DNase_footprint.pdf', 
       height = 3.5, width = 6)

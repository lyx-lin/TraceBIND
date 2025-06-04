library(tidyverse)
source("~/codes/footprints/footprint_prediction.R")
library(ggplot2)
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
index = 4000:8000

### 10x ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects")
pbmc_scatac = read.table("PBMC_atac/obsBias.tsv")
pbmc_multiome = read.table("pbmc_multiome/obsBias.tsv")
brain_multi = read.table("human_brain_scmultiome/obsBias.tsv")
lymph_node = read.table("lymph_node/obsBias.tsv")

### gene silencing  ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/degron")
degron_ctcf_control = read.csv('CTCF/control/mtDNA/obsBias.tsv', sep = '\t')
degron_ctcf_treatment = read.csv('CTCF/treatment/mtDNA/obsBias.tsv', sep = '\t')
degron_polr2a_control = read.csv('POLR2A/control/mtDNA/obsBias.tsv', sep = '\t')
degron_polr2a_treatment = read.csv('POLR2A/treatment/mtDNA/obsBias.tsv', sep = '\t')

### snATAC ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/")
obsbias_SAMN18736215 = read.table('SAMN18736215/mtDNA/obsBias.tsv')
obsbias_SAMN18736216 = read.table('SAMN18736216/mtDNA/obsBias.tsv')
obsbias_SAMN27505541 = read.table('SAMN27505541/mtDNA/obsBias.tsv')
obsbias_SAMN27505542 = read.table('SAMN27505542/mtDNA/obsBias.tsv')
obsbias_SAMN27505543 = read.table('SAMN27505543/mtDNA/obsBias.tsv')
obsbias_SAMN27505544 = read.table('SAMN27505544/mtDNA/obsBias.tsv')

obsbias_Control_1 = read.table('Control_1/mtDNA/obsBias.tsv')
obsbias_Control_2 = read.table('Control_2/mtDNA/obsBias.tsv')
obsbias_Control_3 = read.table('Control_3/mtDNA/obsBias.tsv')
obsbias_Control_4 = read.table('Control_4/mtDNA/obsBias.tsv')
obsbias_Control_5 = read.table('Control_5/mtDNA/obsBias.tsv')
obsbias_Control_6 = read.table('Control_6/mtDNA/obsBias.tsv')

obsbias_DN_1 = read.table('DN_1/mtDNA/obsBias.tsv')
obsbias_DN_2 = read.table('DN_2/mtDNA/obsBias.tsv')
obsbias_DN_3 = read.table('DN_3/mtDNA/obsBias.tsv')
obsbias_DN_4 = read.table('DN_4/mtDNA/obsBias.tsv')
obsbias_DN_5 = read.table('DN_5/mtDNA/obsBias.tsv')
obsbias_DN_6 = read.table('DN_6/mtDNA/obsBias.tsv')

obsbias_CKD_1 = read.table('CKD_1/mtDNA/obsBias.tsv')
obsbias_CKD_2 = read.table('CKD_2/mtDNA/obsBias.tsv')
obsbias_CKD_3 = read.table('CKD_3/mtDNA/obsBias.tsv')
obsbias_CKD_4 = read.table('CKD_4/mtDNA/obsBias.tsv')
obsbias_CKD_5 = read.table('CKD_5/mtDNA/obsBias.tsv')

### Multiome ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")

obsbias_090922Nx = read.table('090922Nx/mtDNA/obsBias.tsv')
obsbias_091422Nx = read.table('091422Nx/mtDNA/obsBias.tsv')
obsbias_1_27Nx = read.table('1-27Nx/mtDNA/obsBias.tsv')
obsbias_2_15Nx = read.table('2-15Nx/mtDNA/obsBias.tsv')
obsbias_A2 = read.table('A2/mtDNA/obsBias.tsv')
obsbias_AIIM164 = read.table('AIIM164/mtDNA/obsBias.tsv')
obsbias_AIL5160 = read.table('AIL5160/mtDNA/obsBias.tsv')
obsbias_AJDL105 = read.table('AJDL105/mtDNA/obsBias.tsv')
obsbias_AJDV174 = read.table('AJDV174/mtDNA/obsBias.tsv')

### heatmap ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
obsbias = data.frame('polr2a_control' = degron_polr2a_control$obs_bias[index], 
                     # 'polr2a_treatment' = degron_polr2a_treatment$obs_bias[index], 
                     'ctcf_control' = degron_ctcf_control$obs_bias[index], 
                     # 'ctcf_treatment' = degron_ctcf_treatment$obs_bias[index], 
                     # "SAMN27505544" = obsbias_SAMN27505544$obs_bias[index], 
                     # "SAMN18736215" = obsbias_SAMN18736215$obs_bias[index], 
                     "SAMN18736216" = obsbias_SAMN18736216$obs_bias[index], 
                     # "SAMN27505541" = obsbias_SAMN27505541$obs_bias[index], 
                     # "SAMN27505543" = obsbias_SAMN27505543$obs_bias[index], 
                     # '1_27Nx' = obsbias_1_27Nx$obs_bias[index], 
                     # '2_15Nx' = obsbias_2_15Nx$obs_bias[index], 
                     # '091422Nx' = obsbias_091422Nx$obs_bias[index], 
                     # '090922Nx' = obsbias_090922Nx$obs_bias[index],
                     'AIL5160' = obsbias_AIL5160$obs_bias[index], 
                     'AJDV174' = obsbias_AJDV174$obs_bias[index], 
                     'AIIM164' = obsbias_AIIM164$obs_bias[index], 
                     # 'brain_multi' = brain_multi$obs_bias[index], 
                     # "Control_1" = obsbias_Control_1$obs_bias[index], 
                     # "Control_2" = obsbias_Control_2$obs_bias[index], 
                     # "Control_3" = obsbias_Control_3$obs_bias[index], 
                     "DN_1" = obsbias_DN_1$obs_bias[index], 
                     "CKD_3" = obsbias_CKD_3$obs_bias[index], 
                     'AJDL105' = obsbias_AJDL105$obs_bias[index], 
                     'pbmc_multiome' = pbmc_multiome$obs_bias[index], 
                     # "CKD_4" = obsbias_CKD_4$obs_bias[index], 
                     # "CKD_5" = obsbias_CKD_5$obs_bias[index], 
                     # "Control_4" = obsbias_Control_4$obs_bias[index], 
                     "Control_5" = obsbias_Control_5$obs_bias[index], 
                     "Control_6" = obsbias_Control_6$obs_bias[index], 
                     "DN_5" = obsbias_DN_5$obs_bias[index], 
                     "CKD_2" = obsbias_CKD_2$obs_bias[index], 
                     # "DN_6" = obsbias_DN_6$obs_bias[index], 
                     # "DN_2" = obsbias_DN_2$obs_bias[index], 
                     # "DN_3" = obsbias_DN_3$obs_bias[index], 
                     # "CKD_1" = obsbias_CKD_1$obs_bias[index], 
                     # 'A2' = obsbias_A2$obs_bias[index],
                     "SAMN27505541" = obsbias_SAMN27505541$obs_bias[index], 
                     'pbmc_scatac' = pbmc_scatac$obs_bias[index], 
                     'lymph_node' = lymph_node$obs_bias[index]
                    )
cor_plot = cor
cor_data = cor(obsbias) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample1') %>%
  pivot_longer(cols = -sample1, 
               names_to = 'sample2', 
               values_to = 'correlation') %>%
  arrange(sample1) %>%
  group_by(sample1) %>%
  filter(row_number() >= which(sample1 == sample2))

cor_data$sample1 = factor(cor_data$sample1, 
                          levels = colnames(obsbias), 
                          labels = gsub("X", "", colnames(obsbias)))
cor_data$sample2 = factor(cor_data$sample2, 
                          levels = rev(colnames(obsbias)), 
                          labels = gsub("X", "", rev(colnames(obsbias))))
# cor_data$correlation[cor_data$correlation <= 0.8] = 0.81
cor_heatmap = ggplot(cor_data, aes(sample1, sample2)) +
  theme_bw() +
  xlab('sample 1') +
  ylab('sample 2') +
  geom_tile(aes(fill = correlation), color='white') +
  scale_fill_gradient(low = 'white', 
                      high = 'darkblue', 
                      space = 'Lab', 
                      limits = c(0.85, 1)
  ) +
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    # axis.text.y=element_text(size = 6),
    # axis.text.x=element_text(angle=90, size = 6),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(), 
    # legend.position='none'
  )
cor_heatmap
ggsave("../cor_heatmap_obs_bias.pdf", 
       cor_heatmap, height = 5, 
       width = 6)

p1 = GGally::ggpairs(obsbias)
ggsave("../cor_scatter_obs_bias.pdf", 
       p1, height = 14, 
       width = 14)
ggsave("../cor_scatter_obs_bias.png", 
       p1, height = 14, 
       width = 14)

### density plot -----
## same study 
data = data.frame('sample1' = obsbias_AJDL105$obs_bias, 
                  'sample2' = obsbias_DN_1$obs_bias)
ggplot(data = data[index, ], 
       aes(x = sample1, y = sample2)) +
  ylab("sample 1") + 
  xlab("sample 2") + 
  geom_hdr(aes(fill = after_stat(probs)), 
           color = "black", alpha = 0.8, 
           probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", 
           label.x = -0.5, 
           digits = 3, 
           label.y = 12) +
  xlim(c(-0.5, 12)) +
  ylim(c(-0.5, 12)) +
  # ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "Blues")) +
  theme_Publication()
ggsave("../obs_bias_same_study.pdf", 
       height = 5, 
       width = 5)

### same lab 
data = data.frame('sample1' = obsbias_AIIM164$obs_bias, 
                  'sample2' = obsbias_CKD_3$obs_bias)
ggplot(data = data[index, ], 
       aes(x = sample1, y = sample2)) +
  ylab("sample 1") + 
  xlab("sample 2") + 
  geom_hdr(aes(fill = after_stat(probs)), 
           color = "black", alpha = 0.8, 
           probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", 
           label.x = -0.5, 
           digits = 3, 
           label.y = 18) +
  # xlim(c(-0.5, 12)) +
  # ylim(c(-0.5, 12)) +
  # ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "Blues")) +
  theme_Publication()
ggsave("../obs_bias_same_lab.pdf", 
       height = 5, 
       width = 5)

### different labs
data = data.frame('sample1' = obsbias_Control_6$obs_bias, 
                  'sample2' = lymph_node$obs_bias)
ggplot(data = data[index, ], 
       aes(x = sample1, y = sample2)) +
  ylab("sample 1") + 
  xlab("sample 2") + 
  geom_hdr(aes(fill = after_stat(probs)), 
           color = "black", alpha = 0.8, 
           probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", 
           label.x = -0.5, 
           digits = 3, 
           label.y = 15) +
  xlim(c(-0.5, 12)) +
  ylim(c(-0.5, 15)) +
  # ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "Blues")) +
  theme_Publication()
ggsave("../obs_bias_different_lab.pdf", 
       height = 5, 
       width = 5)

### downsample ----
data = data.frame('sample1' = obsbias_Control_6$obs_bias, 
                  'sample2' = obsbias_A2$obs_bias)
ggplot(data = data[index, ], 
       aes(x = sample1, y = sample2)) +
  ylab("sample 1") + 
  xlab("sample 1") + 
  geom_hdr(aes(fill = after_stat(probs)), 
           color = "black", alpha = 0.8, 
           probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", 
           label.x = -0.5, 
           digits = 3, 
           label.y = 15) +
  # xlim(c(-0.5, 12)) +
  # ylim(c(-0.5, 12)) +
  # ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "Blues")) +
  theme_Publication()

ggsave("../A2_Control_6.pdf", 
       height = 5, 
       width = 5)
set.seed(42)
downsampled_data = data.frame(matrix(data = 0, 50, 1))
downsampled_data = downsampled_data[-1]

samples = list('polr2a_control' = degron_polr2a_control, 
                     'polr2a_treatment' = degron_polr2a_treatment, 
                     # 'ctcf_control' = degron_ctcf_control, 
                     # 'ctcf_treatment' = degron_ctcf_treatment, 
                     # "SAMN27505544" = obsbias_SAMN27505544, 
                     # "SAMN18736215" = obsbias_SAMN18736215, 
                     "SAMN18736216" = obsbias_SAMN18736216, 
                     # "SAMN27505541" = obsbias_SAMN27505541, 
                     # "SAMN27505543" = obsbias_SAMN27505543, 
                     # '1_27Nx' = obsbias_1_27Nx, 
                     # '2_15Nx' = obsbias_2_15Nx, 
                     # '091422Nx' = obsbias_091422Nx, 
                     # '090922Nx' = obsbias_090922Nx,
                     'AIL5160' = obsbias_AIL5160, 
                     'AJDV174' = obsbias_AJDV174, 
                     'AIIM164' = obsbias_AIIM164, 
                     # 'brain_multi' = brain_multi, 
                     # "Control_1" = obsbias_Control_1, 
                     # "Control_2" = obsbias_Control_2, 
                     # "Control_3" = obsbias_Control_3, 
                     "DN_1" = obsbias_DN_1, 
                     "CKD_3" = obsbias_CKD_3, 
                     'AJDL105' = obsbias_AJDL105, 
                     'pbmc_multiome' = pbmc_multiome, 
                     # "CKD_4" = obsbias_CKD_4, 
                     # "CKD_5" = obsbias_CKD_5, 
                     # "Control_4" = obsbias_Control_4, 
                     "Control_5" = obsbias_Control_5, 
                     "Control_6" = obsbias_Control_6, 
                     "DN_5" = obsbias_DN_5, 
                     "CKD_2" = obsbias_CKD_2, 
                     # "DN_6" = obsbias_DN_6, 
                     # "DN_2" = obsbias_DN_2, 
                     # "DN_3" = obsbias_DN_3, 
                     # "CKD_1" = obsbias_CKD_1, 
                     # 'dabtram_day10' = obsbias_dabtram_day10, 
                     # 'cis_day10' = obsbias_cis_day10, 
                     # 'A2' = obsbias_A2,
                     'dabtram_week5' = obsbias_dabtram_week5,
                     'cis_week5' = obsbias_cis_week5, 
                     'day0' = obsbias_day0, 
                     "SAMN27505541" = obsbias_SAMN27505541, 
                     'pbmc_scatac' = pbmc_scatac, 
                     'lymph_node' = lymph_node
)

sample_pairs <- combn(seq_along(samples), 2, simplify = FALSE)

# Function to compute downsampled correlations for one pair
process_pair <- function(pair) {
  i <- pair[1]
  j <- pair[2]
  data_1 <- samples[[i]]
  data_2 <- samples[[j]]
  
  if ('context' %in% colnames(data_1)) data_1$context <- NULL
  if ('context' %in% colnames(data_2)) data_2$context <- NULL
  
  ratio <- min(sum(data_1$insertions)/sum(data_2$insertions),
               sum(data_2$insertions)/sum(data_1$insertions))
  ratio = max(0.3, ratio)
  downsampled_cor <- numeric(50)
  
  if (sum(data_1$insertions) >= sum(data_2$insertions)) {
    for (k in 1:50) {
      set.seed(k)
      data_1$downsampled_insertions <- rbinom(n = nrow(data_1), 
                                              size = data_1$insertions, 
                                              prob = ratio)
      data_1$obs_bias_downsampled <- mapply(function(pos, val) {
        window <- data_1[abs(data_1$positions - pos) <= 50, "downsampled_insertions"]
        val / mean(window)
      }, data_1$positions, data_1[, 5])
      
      downsampled_cor[k] <- cor(data_2$obs_bias[index], 
                                data_1$obs_bias_downsampled[index])
      
    }
  } else {
    for (k in 1:50) {
      set.seed(k)
      data_2$downsampled_insertions <- rbinom(n = nrow(data_2), 
                                              size = data_2$insertions, 
                                              prob = ratio)
      data_2$obs_bias_downsampled <- mapply(function(pos, val) {
        window <- data_2[abs(data_2$positions - pos) <= 50, "downsampled_insertions"]
        val / mean(window)
      }, data_2$positions, data_2[, 5])
      
      downsampled_cor[k] <- cor(data_1$obs_bias[index], 
                                data_2$obs_bias_downsampled[index], 
                                use = 'complete.obs')
    }
  }
  
  return(setNames(list(downsampled_cor), paste0(i, "_", j)))
}

# Run in parallel
results <- pbmcapply::pbmclapply(sample_pairs, process_pair, 
                                 mc.cores = 2)

# Combine results into a named list
downsampled_data <- do.call(c, results)
write.table(downsampled_data, 
            '../downsampled.txt')

for (i in 1:(length(samples) - 1)){
  print(i)
  for (j in (i+1):length(samples)){
    data_1 <- samples[[i]]
    data_2 <- samples[[j]]
    downsampled = downsampled_data[paste0(i, '_', j)]
    data1 = data.frame('y1' = data.frame(downsampled)[, 1], 
                       'x1' = 'downsampled')
    obs_cor = cor(data_1$obs_bias, data_2$obs_bias, 
                  use = 'complete.obs')
    p2 = ggplot() + 
      geom_boxplot(data = data1, aes(x = x1, y = round(y1, 2))) + 
      geom_point(data = data.frame(x1 = 'downsampled', 
                                   y1 = round(obs_cor, 2)),
                 aes(x = x1, y = y1),
                 color = 'red', size = 2) + 
      xlab("") + 
      ylab("") + 
      scale_y_continuous(limits = c(0.8, 1), 
                         breaks = seq(0.8, 1.0, 0.05))
    p1[j, i] = p2
  }
}
ggsave("../downsample_plot.png", p1, height = 14, width = 14)

### GC content ----
obsbias_Control_1$GC_content = sapply(1:dim(obsbias_Control_1)[1], function(x) 
  sum(strsplit(obsbias_Control_1$context[x], "")[[1]] %in% c("C", "G"))/101)

obsbias = data.frame('polr2a_control' = degron_polr2a_control$obs_bias, 
                     # 'polr2a_treatment' = degron_polr2a_treatment$obs_bias, 
                     'ctcf_control' = degron_ctcf_control$obs_bias, 
                     # 'ctcf_treatment' = degron_ctcf_treatment$obs_bias, 
                     # "SAMN27505544" = obsbias_SAMN27505544$obs_bias, 
                     # "SAMN18736215" = obsbias_SAMN18736215$obs_bias, 
                     "SAMN18736216" = obsbias_SAMN18736216$obs_bias, 
                     # "SAMN27505541" = obsbias_SAMN27505541$obs_bias, 
                     # "SAMN27505543" = obsbias_SAMN27505543$obs_bias, 
                     # '1_27Nx' = obsbias_1_27Nx$obs_bias, 
                     # '2_15Nx' = obsbias_2_15Nx$obs_bias, 
                     # '091422Nx' = obsbias_091422Nx$obs_bias, 
                     # '090922Nx' = obsbias_090922Nx$obs_bias,
                     'AIL5160' = obsbias_AIL5160$obs_bias, 
                     'AJDV174' = obsbias_AJDV174$obs_bias, 
                     'AIIM164' = obsbias_AIIM164$obs_bias, 
                     # 'brain_multi' = brain_multi$obs_bias, 
                     # "Control_1" = obsbias_Control_1$obs_bias, 
                     # "Control_2" = obsbias_Control_2$obs_bias, 
                     # "Control_3" = obsbias_Control_3$obs_bias, 
                     "DN_1" = obsbias_DN_1$obs_bias, 
                     "CKD_3" = obsbias_CKD_3$obs_bias, 
                     'AJDL105' = obsbias_AJDL105$obs_bias, 
                     'pbmc_multiome' = pbmc_multiome$obs_bias, 
                     # "CKD_4" = obsbias_CKD_4$obs_bias, 
                     # "CKD_5" = obsbias_CKD_5$obs_bias, 
                     # "Control_4" = obsbias_Control_4$obs_bias, 
                     "Control_5" = obsbias_Control_5$obs_bias, 
                     "Control_6" = obsbias_Control_6$obs_bias, 
                     "DN_5" = obsbias_DN_5$obs_bias, 
                     "CKD_2" = obsbias_CKD_2$obs_bias, 
                     # "DN_6" = obsbias_DN_6$obs_bias, 
                     # "DN_2" = obsbias_DN_2$obs_bias, 
                     # "DN_3" = obsbias_DN_3$obs_bias, 
                     # "CKD_1" = obsbias_CKD_1$obs_bias, 
                     # 'dabtram_day10' = obsbias_dabtram_day10$obs_bias, 
                     # 'cis_day10' = obsbias_cis_day10$obs_bias, 
                     # 'A2' = obsbias_A2$obs_bias,
                     'dabtram_week5' = obsbias_dabtram_week5$obs_bias,
                     'cis_week5' = obsbias_cis_week5$obs_bias, 
                     'day0' = obsbias_day0$obs_bias, 
                     "SAMN27505541" = obsbias_SAMN27505541$obs_bias, 
                     'pbmc_scatac' = pbmc_scatac$obs_bias, 
                     'lymph_node' = lymph_node$obs_bias
)
obsbias$GC_content = obsbias_Control_1$GC_content
obsbias = obsbias %>% 
  pivot_longer(!GC_content,  
               names_to = 'dataset',         
               values_to = 'insertions')
obsbias = obsbias %>%
  group_by(GC_content, dataset) %>%
  summarise(insertions = mean(insertions))

obsbias_mean <- obsbias %>%
  group_by(GC_content) %>%
  summarise(
    mean = mean(insertions),
    sd = sd(insertions)
  )
summerNight = c("#78B7C5", "#F49B7C", "#fbdf72",
                "#E6C2DC","#9cdff0","#fbbf45","#bed678")

ggplot(obsbias_mean) + 
  # geom_point(aes(x = GC_content, y = bias_diff)) + 
  geom_ribbon(aes(x = GC_content, 
                  ymin = mean - 1.96*sd, 
                  ymax = mean + 1.96*sd, 
                  fill = "1.96*Standard Deviation"), 
              alpha = 0.5) +
  ylab("observed bias") + 
  geom_smooth(aes(x = GC_content, y = mean, fill = "Mean"), 
              span = 0.1, se = F, color = '#064ACB') + 
  scale_fill_manual(name = '',  
                    values = c("1.96*Standard Deviation" = "#FFB200", 
                               'Mean' = '#064ACB')) +
  theme_classic() + 
  theme(legend.position=c(0.7, 0.9)) + 
  theme_Publication()

ggsave("../gc_contents_insertions1.pdf", 
       height = 3, width = 3)

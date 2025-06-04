source("~/codes/footprints/footprint_prediction.R")
library(tidyverse)
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
theme_Publication<- function(base_size=12, base_family="sans") {
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
## multiome -----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")

samples_multiome_names = c("multiome_kidney_mtDNA/1-27Nx", 
                           "multiome_kidney_mtDNA/2-15Nx",
                           "multiome_kidney_mtDNA/090922Nx", 
                           "multiome_kidney_mtDNA/091422Nx", 
                           "multiome_kidney_mtDNA/A2", 
                           "multiome_kidney_mtDNA/AIIM164", 
                           "multiome_kidney_mtDNA/AIL5160", 
                           "multiome_kidney_mtDNA/AJDL105", 
                           "multiome_kidney_mtDNA/AJDV174"
)
count_090922Nx = readRDS('090922Nx/mtDNA/counts.rds')
count_091422Nx = readRDS('091422Nx/mtDNA/counts.rds')
count_1_27Nx = readRDS('1-27Nx/mtDNA/counts.rds')
count_2_15Nx = readRDS('2-15Nx/mtDNA/counts.rds')
count_A2 = readRDS('A2/mtDNA/counts.rds')
count_AIIM164 = readRDS('AIIM164/mtDNA/counts.rds')
count_AIL5160 = readRDS('AIL5160/mtDNA/counts.rds')
count_AJDL105 = readRDS('AJDL105/mtDNA/counts.rds')
count_AJDV174 = readRDS('AJDV174/mtDNA/counts.rds')

samples_multiome = list(count_1_27Nx, 
                        count_2_15Nx, 
                        count_090922Nx,
                        count_091422Nx, 
                        count_A2, 
                        count_AIIM164, 
                        count_AIL5160,
                        count_AJDL105, 
                        count_AJDV174
)

celltype_kidney = c("DCT1", "DCT2", "PC", "PCT", "PEC", "ICA", 
                    "ICB", "PODO", "PST", 
                    "PT_VCAM1", "TAL1","TAL2")
celltype_others = c("BCELL", 
                    "ENDO", 
                    "FIB_VSMC_MC", 
                    "MONO", 
                    "TCELL")

for (j in 1:length(samples_multiome)){
  counts = samples_multiome[[j]]
  celltype_bias = data.frame(matrix(data = 0, 16000, 1))
  celltype_bias = celltype_bias[-1]
  
  insertions_kidney = NULL
  for (i in 1:16){
    insertions_kidney = c(insertions_kidney, get_insertions(countData = counts, 
                                                            regionInd = i, 
                                                            groupIDs = celltype_kidney, 
                                                            width = 1000))
  }
  obsbias = data.frame('positions' = 501:16500, 
                       'insertions' = as.numeric(insertions_kidney))
  obsbias$insertions = as.numeric(obsbias$insertions)
  obsbias$positions = as.numeric(obsbias$positions)
  
  obsbias_kidney = apply(obsbias, 1, function(x){
    as.numeric(x[2])/as.numeric(mean(obsbias[as.numeric(obsbias$positions) >= as.numeric(x[1]) - 50 & 
                                               as.numeric(obsbias$positions) <= as.numeric(x[1]) + 50, "insertions"]))
  })
  celltype_bias['kidney_cells'] = obsbias_kidney
  insertions_others = NULL
  for (i in 1:16){
    insertions_others = c(insertions_others, get_insertions(countData = counts, 
                                                            regionInd = i, 
                                                            groupIDs = celltype_others, 
                                                            width = 1000))
  }
  obsbias = data.frame('positions' = 501:16500, 
                       'insertions' = as.numeric(insertions_others))
  obsbias$insertions = as.numeric(obsbias$insertions)
  obsbias$positions = as.numeric(obsbias$positions)
  
  obsbias_other = apply(obsbias, 1, function(x){
    as.numeric(x[2])/as.numeric(mean(obsbias[as.numeric(obsbias$positions) >= as.numeric(x[1]) - 50 & 
                                               as.numeric(obsbias$positions) <= as.numeric(x[1]) + 50, "insertions"]))
  })
  celltype_bias['other_cells'] = obsbias_other
  cor_num = round(cor(celltype_bias$kidney_cells, 
                  celltype_bias$other_cells), 2)
  p1 <- ggplot(celltype_bias, 
               aes(x = other_cells, y = kidney_cells)) +
    ylab("kidney cells") + 
    xlab("other cells") + 
    geom_point() + 
    annotate("text", x = 4, y = 12, 
             label = paste0("R = ", cor_num, ', p < 2.2e-16'), 
             color = 'red', size = 6) + 
    theme_Publication()
  ggsave(paste0("../", samples_multiome_names[j], 
                "/", strsplit(samples_multiome_names[j], '/')[[1]][2], "_obs_bias_kidney_v.s._others_scatter.pdf"), 
         p1, height = 5, width = 5)
}

## snATAC ----
samples_names_snatac = c("snATAC_mtDNA/Control_1", 
                         "snATAC_mtDNA/Control_2", 
                         "snATAC_mtDNA/Control_3", 
                         "snATAC_mtDNA/Control_4", 
                         "snATAC_mtDNA/Control_5",
                         "snATAC_mtDNA/Control_6",
                         "snATAC_mtDNA/DN_1", 
                         "snATAC_mtDNA/DN_2", 
                         "snATAC_mtDNA/DN_3", 
                         "snATAC_mtDNA/DN_4", 
                         "snATAC_mtDNA/DN_5", 
                         "snATAC_mtDNA/CKD_1", 
                         "snATAC_mtDNA/CKD_2", 
                         "snATAC_mtDNA/CKD_3", 
                         "snATAC_mtDNA/CKD_4", 
                         "snATAC_mtDNA/CKD_5", 
                         "snATAC_mtDNA/SAMN18736215", 
                         "snATAC_mtDNA/SAMN18736216", 
                         "snATAC_mtDNA/SAMN27505541", 
                         "snATAC_mtDNA/SAMN27505542", 
                         "snATAC_mtDNA/SAMN27505543", 
                         "snATAC_mtDNA/SAMN27505544")

setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/")

count_Control_1 = readRDS('Control_1/mtDNA/counts.rds')
count_Control_2 = readRDS('Control_2/mtDNA/counts.rds')
count_Control_3 = readRDS('Control_3/mtDNA/counts.rds')
count_Control_4 = readRDS('Control_4/mtDNA/counts.rds')
count_Control_5 = readRDS('Control_5/mtDNA/counts.rds')
count_Control_6 = readRDS('Control_6/mtDNA/counts.rds')

count_DN_1 = readRDS('DN_1/mtDNA/counts.rds')
count_DN_2 = readRDS('DN_2/mtDNA/counts.rds')
count_DN_3 = readRDS('DN_3/mtDNA/counts.rds')
count_DN_4 = readRDS('DN_4/mtDNA/counts.rds')
count_DN_5 = readRDS('DN_5/mtDNA/counts.rds')

count_CKD_1 = readRDS('CKD_1/mtDNA/counts.rds')
count_CKD_2 = readRDS('CKD_2/mtDNA/counts.rds')
count_CKD_3 = readRDS('CKD_3/mtDNA/counts.rds')
count_CKD_4 = readRDS('CKD_4/mtDNA/counts.rds')
count_CKD_5 = readRDS('CKD_5/mtDNA/counts.rds')

count_SAMN18736215 = readRDS('SAMN18736215/mtDNA/counts.rds')
count_SAMN18736216 = readRDS('SAMN18736216/mtDNA/counts.rds')
count_SAMN27505541 = readRDS('SAMN27505541/mtDNA/counts.rds')
count_SAMN27505542 = readRDS('SAMN27505542/mtDNA/counts.rds')
count_SAMN27505543 = readRDS('SAMN27505543/mtDNA/counts.rds')
count_SAMN27505544 = readRDS('SAMN27505544/mtDNA/counts.rds')

samples_snATAC = list(count_Control_1, 
                      count_Control_2, 
                      count_Control_3,
                      count_Control_4, 
                      count_Control_5, 
                      count_Control_6,
                      count_DN_1, 
                      count_DN_2, 
                      count_DN_3,
                      count_DN_4, 
                      count_DN_5, 
                      count_CKD_1, 
                      count_CKD_2, 
                      count_CKD_3,
                      count_CKD_4, 
                      count_CKD_5, 
                      count_SAMN18736215, 
                      count_SAMN18736216, 
                      count_SAMN27505541, 
                      count_SAMN27505542, 
                      count_SAMN27505543, 
                      count_SAMN27505544)

unique(count_CKD_1[[1]]$group)
celltype_kidney = c("DCT1", "DCT2_PC",
                    "PCT", "PST", "TAL",
                    "PT_PROM1", "PT_VCAM1")
celltype_others = c("BCELL", "ENDO", "MONO", "TCELL", 
                    "FIB_VSMC_MC", "PEC", "PODO", 
                    "ICA", "ICB")
for (j in 1:length(samples_snATAC)){
  counts = samples_snATAC[[j]]
  celltype_bias = data.frame(matrix(data = 0, 16000, 1))
  celltype_bias = celltype_bias[-1]
  
  insertions_kidney = NULL
  for (i in 1:16){
    insertions_kidney = c(insertions_kidney, get_insertions(countData = counts, 
                                                            regionInd = i, 
                                                            groupIDs = celltype_kidney, 
                                                            width = 1000))
  }
  obsbias = data.frame('positions' = 501:16500, 
                       'insertions' = as.numeric(insertions_kidney))
  obsbias$insertions = as.numeric(obsbias$insertions)
  obsbias$positions = as.numeric(obsbias$positions)
  
  obsbias_kidney = apply(obsbias, 1, function(x){
    as.numeric(x[2])/as.numeric(mean(obsbias[as.numeric(obsbias$positions) >= as.numeric(x[1]) - 50 & 
                                               as.numeric(obsbias$positions) <= as.numeric(x[1]) + 50, "insertions"]))
  })
  celltype_bias['kidney_cells'] = obsbias_kidney
  insertions_others = NULL
  for (i in 1:16){
    insertions_others = c(insertions_others, get_insertions(countData = counts, 
                                                            regionInd = i, 
                                                            groupIDs = celltype_others, 
                                                            width = 1000))
  }
  obsbias = data.frame('positions' = 501:16500, 
                       'insertions' = as.numeric(insertions_others))
  obsbias$insertions = as.numeric(obsbias$insertions)
  obsbias$positions = as.numeric(obsbias$positions)
  
  obsbias_other = apply(obsbias, 1, function(x){
    as.numeric(x[2])/as.numeric(mean(obsbias[as.numeric(obsbias$positions) >= as.numeric(x[1]) - 50 & 
                                               as.numeric(obsbias$positions) <= as.numeric(x[1]) + 50, "insertions"]))
  })
  celltype_bias['other_cells'] = obsbias_other
  
  p1 <- ggplot(celltype_bias[index, ], 
               aes(x = other_cells, y = kidney_cells)) +
    ylab("kidney cells") + 
    xlab("Other cells") + 
    geom_point() + 
    # geom_hdr(aes(fill = after_stat(probs)), 
    #          color = "black", alpha = 0.8, 
    #          probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
    annotate("text", x = 4, y = 12, 
             label = paste0("R = ", cor_num, ', p < 2.2e-16'), 
             color = 'red', size = 6) +
    theme_Publication()
  cor(celltype_bias[, 'other_cells'], 
      celltype_bias[, 'kidney_cells'])
  ggsave(paste0("../", samples_names_snatac[j], "/", strsplit(samples_names_snatac[j], '/')[[1]][2], "_obs_bias_kidney_v.s._others_scatter.pdf"), 
         p1, height = 5, width = 5)
}


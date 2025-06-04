args <- commandArgs(TRUE)
x = as.numeric(args[1])
k = as.numeric(args[2])

print(x)

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reticulate)
library(tensorflow)
library(Biostrings)
library(GenomicRanges)

projectNames = c('CTCF/treatment', 
                 'CTCF/control')
projectName = projectNames[k]
print(projectName)
setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/degron/", projectName))

source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("../peaks_region.txt")
colnames(regionsBed) = c('chr', 'start', 'end')

runs = (100*x-99):min(dim(regionsBed)[1], 100*x)

bias = read.table('freezed_finetuned_pred_bias.txt')
counts = readRDS("counts.rds")

for (i in runs){
  print(i)
  bardata = get_bardata(counts = counts, 
                        regionsBed = regionsBed, 
                        bias = bias[i, 
                                    1:length(regionsBed$start[i]:regionsBed$end[i])], 
                        i = i, 
                        groupIDs = '1')

  Tn5_plot = plot_insertions(bardata$Tn5Insertion, 
                             positions = bardata$position, 
                             chr = regionsBed$chr[i])
  
  footprinting_results = NB_footprintings(Tn5Insertion = bardata$Tn5Insertion, 
                                          pred_bias = bardata$pred_bias, 
                                          positions = bardata$position, 
                                          p.adjust.method = 'BH', 
                                          nCores = 2)
  
  p_value_matrix = footprinting_results[['pval']]
  effect_size_matrix = footprinting_results[['effect_size']]
  
  if(!dir.exists(paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/", projectName))){
    system(paste("mkdir -p", paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/", projectName)))
  }
  
  write.table(effect_size_matrix, 
              paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/", projectName, '/', 
                     i, "_effect_size_matrix.txt"))
  write.table(p_value_matrix, 
              paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/", projectName, '/', 
                     i, "_p_value_matrix.txt"))
  
}


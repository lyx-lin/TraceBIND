args <- commandArgs(TRUE)
x = as.numeric(args[[1]])
print(x)
runs = (70*x-69):(70*x)

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reticulate)
library(tensorflow)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

projectName = "kidney_benchmarking/NR3C1/human_kidney/Control6"

if(!dir.exists(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/", projectName))){
  system(paste("mkdir -p", paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/", projectName)))
}

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/", projectName))

source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("peaks_region.txt")
bias_Control6 = read.table('freezed_finetuned_pred_bias.txt')
counts = readRDS("../aggregated_counts.rds")
names(counts) = as.character(1:length(counts))

for (i in runs){
  print(i)
  pred_bias = bias_Control6[i, 1:length(regionsBed$start[i]:regionsBed$end[i])]
  bardata = get_bardata(counts = counts,  
                        regionsBed = regionsBed, 
                        bias = pred_bias, 
                        i = i, 
                        groupIDs = unique(counts[[i]]$group))
  
  Tn5_plot = plot_insertions(bardata$Tn5Insertion, 
                             bardata$position, 
                             regionsBed$chr[i])
  
  footprinting_results = NB_footprintings(bardata$Tn5Insertion, 
                                          bardata$pred_bias, 
                                          bardata$position, 
                                          p.adjust.method = 'BH', 
                                          nCores = 2)
  
  p_value_matrix = footprinting_results[['pval']]
  effect_size_matrix = footprinting_results[['effect_size']]
  
  source("~/codes/footprints/footprint_prediction.R")
  binding_sites = binding_sites(p_value_matrix, 
                                effect_size_matrix, 
                                width_threshold = 0, 
                                p_threshold = 0.1)
  if (dim(binding_sites)[1] > 0){
    binding_sites$coverage = unlist(sapply(1:dim(binding_sites)[1], function(index) {
      positions_1 <- (binding_sites$position[index] - binding_sites$width[index] / 2 - 1):
        (binding_sites$position[index] - 1 - binding_sites$width[index] * 3 / 2)
      positions_2 <- (binding_sites$position[index] + binding_sites$width[index] / 2 + 1):
        (binding_sites$position[index] + 1 + binding_sites$width[index] * 3 / 2)
      
      mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                          mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
      mean_coverage
    }))
  }
  
  if(!dir.exists(paste0("~/nzhanglab/project/linyx/footprints/results/data/NR3C1/effect_size/nb_bs_BH"))){
    system(paste("mkdir -p", paste0("~/nzhanglab/project/linyx/footprints/results/data/NR3C1/effect_size/nb_bs_BH")))
  }
  
  write.table(binding_sites, 
              paste0("~/nzhanglab/project/linyx/footprints/results/data/NR3C1/effect_size/nb_bs_BH/",
                     i, "_", "binding_sites_nb.txt"))
  
}

args <- commandArgs(TRUE)
x = as.numeric(args[1])
k = as.numeric(args[2])
print(x)
print(k)

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reticulate)
library(tensorflow)
library(Biostrings)
library(GenomicRanges)

projectNames = c('age_82_weeks/U21/kidney_multiome', 
                 'age_82_weeks/U21/kidney_defnd', 
                 'age_82_weeks/V22/kidney_defnd')
projectName = projectNames[k]

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName))

source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("../../peaks_region.txt")
colnames(regionsBed) = c('chr', 'start', 'end')

runs = (100*x-99):min(dim(regionsBed)[1], 100*x)

bias = read.table('freezed_finetuned_pred_bias.txt', 
                  nrows = 100, skip = 100*x - 100)

counts = readRDS("counts_subset.rds")
counts = counts[runs]
regionsBed = regionsBed[runs, ]

for (i in runs){
  print(i)
  bardata = get_bardata(counts = counts, 
                        regionsBed = regionsBed, 
                        bias = bias[i - (100*x - 99) + 1, 
                                    1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                        i = i - (100*x - 99) + 1, 
                        groupIDs = '1')

  Tn5_plot = plot_insertions(bardata$Tn5Insertion, 
                             positions = bardata$position, 
                             chr = regionsBed$chr[ i - (100*x - 99) + 1])
  
  footprinting_results = NB_footprintings(Tn5Insertion = bardata$Tn5Insertion, 
                                          pred_bias = bardata$pred_bias, 
                                          positions = bardata$position, 
                                          p.adjust.method = 'BH', 
                                          nCores = 2)
  
  p_value_matrix = footprinting_results[['pval']]
  effect_size_matrix = footprinting_results[['effect_size']]
  
  source("~/codes/footprints/footprint_prediction.R")
  binding_sites = binding_sites(p_value_matrix, 
                                effect_size_matrix,
                                width_threshold = 10, 
                                p_threshold = 0.05)
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
  # binding_sites = binding_sites[binding_sites$effect_size >= 1e-4, ]
  if(!dir.exists(paste0("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/", projectName, "_subset_effect_size_freezed"))){
    system(paste("mkdir -p", paste0("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/", projectName, "_subset_effect_size_freezed")))
  }
  
  write.table(binding_sites, 
              paste0("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND/", projectName,  "_subset_effect_size_freezed/", 
                     i, "_binding_sites_nb.txt"))
  
}

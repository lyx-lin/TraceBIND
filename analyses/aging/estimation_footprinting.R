args <- commandArgs(TRUE)
x = as.numeric(args[1])
k = as.numeric(args[2])
print(x)
# print(k)

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
# library(reticulate)
# library(tensorflow)
library(Biostrings)
library(GenomicRanges)

projectNames = c('age_16_weeks/C03/kidney_multiome', 
                 'age_30_weeks/L12/kidney_multiome',
                 'age_56_weeks/Q17/kidney_multiome', 
                 'age_82_weeks/U21/kidney_multiome')
count_name = c("counts.rds", 
               "counts_subset.rds", 
               "counts_subset.rds", 
               "counts_subset.rds")
projectName = projectNames[k]
count_name = count_names[k]
print(c(projectName, count_name))
setwd(paste0("/~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/consensus/", projectName))
source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("../../../peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')
runs = (100*x-99):min(dim(regionsBed)[1], 100*x)

bias = read.table("freezed_finetuned_pred_bias.txt", 
                  nrows = 100, 
                  skip = 100*x - 100)
counts = readRDS(count_name)
counts = counts[runs]
regionsBed = regionsBed[runs, ]

for (i in runs){
  bardata = get_bardata(counts = counts, 
                        regionsBed = regionsBed, 
                        bias = bias[i - (100*x - 99) + 1,                                         
                                    1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                        i = i - (100*x - 99) + 1, 
                        groupIDs = c("PT", "PT-MT"))
  
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
  
  binding_sites_effect_sizes = binding_sites(p_value_matrix, 
                                             effect_size_matrix, 
                                             p_threshold = 0.05, 
                                             width_threshold = 10)

  binding_sites_effect_sizes['mean_coverage'] = rep(NA, dim(binding_sites_effect_sizes)[1])
  if (dim(binding_sites_effect_sizes)[1] > 0){
    for (index in (1:dim(binding_sites_effect_sizes)[1])){
      positions_1 = (binding_sites_effect_sizes$position[index] - binding_sites_effect_sizes$width[index]/2 - 1):(binding_sites_effect_sizes$position[index] - 1 - binding_sites_effect_sizes$width[index]*3/2)
      positions_2 = (binding_sites_effect_sizes$position[index] + binding_sites_effect_sizes$width[index]/2 + 1):(binding_sites_effect_sizes$position[index] + 1 + binding_sites_effect_sizes$width[index]*3/2)
      mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                          mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
      binding_sites_effect_sizes[index, 'mean_coverage'] = mean_coverage
    }
  }
  if(!dir.exists(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/", projectName, "/effect_size/"))){
    system(paste("mkdir -p", paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/", projectName, "/effect_size/")))
  }
  
  write.table(binding_sites_effect_sizes, 
              paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/consensus/", projectName, "/effect_size/", 
                     i, "_binding_sites_nb.txt"))
  
}

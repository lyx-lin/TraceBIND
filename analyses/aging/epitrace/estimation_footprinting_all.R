args <- commandArgs(TRUE)
x = as.numeric(args[1])
k = as.numeric(args[2])
print(x)
print(k)

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
# library(reticulate)
# library(tensorflow)
library(Biostrings)
library(GenomicRanges)

projectNames = c('age_16_weeks/C03', 
                 'age_30_weeks/L12',
                 'age_56_weeks/Q17', 
                 'age_82_weeks/U21')
setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/"))

source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')

index = min(ceiling(dim(regionsBed)[1]/60000), ceiling(x/600))
bias_16 = read.table(paste0("age_16_weeks/C03/freezed_finetuned_pred_bias_", index, ".txt"), 
                     nrows = 100, 
                     skip = 100*(x - 600*(index-1)) - 100
)
bias_30 = read.table(paste0("age_30_weeks/L12/freezed_finetuned_pred_bias_", index, ".txt"), 
                     nrows = 100, 
                     skip = 100*(x - 600*(index-1)) - 100
)
bias_56 = read.table(paste0("age_56_weeks/Q17/freezed_finetuned_pred_bias_", index, ".txt"), 
                     nrows = 100, 
                     skip = 100*(x - 600*(index-1)) - 100
)
bias_82 = read.table(paste0("age_82_weeks/U21/freezed_finetuned_pred_bias_", index, ".txt"), 
                     nrows = 100, 
                     skip = 100*(x - 600*(index-1)) - 100
)

counts_16 = readRDS("age_16_weeks/C03/counts.rds")
counts_30 = readRDS("age_30_weeks/L12/counts.rds")
counts_56 = readRDS("age_56_weeks/Q17/counts.rds")
counts_82 = readRDS("age_82_weeks/U21/counts.rds")
counts_16 = counts_16[runs]
counts_30 = counts_30[runs]
counts_56 = counts_56[runs]
counts_82 = counts_82[runs]

regionsBed = regionsBed[runs, ]

for (i in runs){
  print(i)
  
  bardata_16 = get_bardata(counts = counts_16, 
                           regionsBed = regionsBed, 
                           bias = bias_16[i - (100*x - 99) + 1, 
                                          1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                           i = i - (100*x - 99) + 1, 
                           groupIDs = c("PT", "PT-MT"))
  bardata_30 = get_bardata(counts = counts_30, 
                           regionsBed = regionsBed, 
                           bias = bias_30[i - (100*x - 99) + 1, 
                                          1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                           i = i - (100*x - 99) + 1, 
                           groupIDs = c("PT", "PT-MT"))
  bardata_56 = get_bardata(counts = counts_56, 
                           regionsBed = regionsBed, 
                           bias = bias_56[i - (100*x - 99) + 1, 
                                          1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                           i = i - (100*x - 99) + 1, 
                           groupIDs = c("PT", "PT-MT"))
  bardata_82 = get_bardata(counts = counts_82, 
                           regionsBed = regionsBed, 
                           bias = bias_82[i - (100*x - 99) + 1,                                         1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
                           i = i - (100*x - 99) + 1, 
                           groupIDs = c("PT", "PT-MT"))
  bardata = bardata_16
  bardata$pred_bias = apply(rbind(bardata_16$pred_bias, 
                                  bardata_30$pred_bias, 
                                  bardata_56$pred_bias, 
                                  bardata_82$pred_bias), 2, 
                            function(x) 
                              weighted.mean(x, c(sum(bardata_16$Tn5Insertion), 
                                                 sum(bardata_30$Tn5Insertion), 
                                                 sum(bardata_56$Tn5Insertion), 
                                                 sum(bardata_82$Tn5Insertion))))
  bardata$Tn5Insertion = bardata_16$Tn5Insertion + bardata_30$Tn5Insertion + 
    bardata_56$Tn5Insertion + bardata_82$Tn5Insertion
  # bardata = get_bardata(counts = counts, 
  #                       regionsBed = regionsBed, 
  #                       bias = bias[i - (100*x - 99) + 1,                                         
  #                                   1:length(regionsBed$start[i - (100*x - 99) + 1]:regionsBed$end[i - (100*x - 99) + 1])], 
  #                       i = i - (100*x - 99) + 1, 
  #                       groupIDs = c("PT", "PT-MT"))
  
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
  if(!dir.exists(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/all_ages/effect_size/"))){
    system(paste("mkdir -p", paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/all_ages/effect_size/")))
  }
  
  write.table(binding_sites_effect_sizes, 
              paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/all_ages/effect_size/", 
                     i, "_binding_sites_nb.txt"))
  
}

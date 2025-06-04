args <- commandArgs(TRUE)
k = as.numeric(args[[1]])
i=11
library(tidyverse)
# reticulate::use_condaenv('base')

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/SAMN27505541"))

source("~/codes/footprints/footprint_prediction.R")
source("../../../../code/utils.R")
source("../../../../code/getBias.R")

regionsBed = data.frame('chr' = 'chrM', 
                        'start' = seq(501, 16500, 1000), 
                        'end' = seq(1500, 16500, 1000))

# get_bias(regionsBed, 
#          referenceGenome = 'hg38', 
#          path = 'mtDNA/finetuned_', 
#          code_path = '../../../../', 
#          model_use = 'Tn5_NN_model_SAMN27505541_finetuned.h5')

# get_bias(regionsBed, 
#          referenceGenome = 'hg38', 
#          path = 'mtDNA/PRINT_', 
#          code_path = '../../../../')

counts = readRDS("mtDNA/counts.rds")
bias_finetuned = read.table("mtDNA/finetuned_pred_bias.txt")
bias_PRINT = read.table("../../PRINT_pred_bias.txt")

pred_bias_finetuned = bias_finetuned[i, ]
pred_bias_PRINT = bias_PRINT[i, ]

ratios = c(seq(0.1, 2, 0.1), 
           seq(2, 16, 1), 
           seq(18, 20, 2))

for (r in ratios){
  set.seed(k)
  ### finetuned ---- 
  bardata_finetuned = get_bardata(counts = counts, 
                                  regionsBed = regionsBed, 
                                  bias = pred_bias_finetuned, 
                                  i = i, 
                                  groupIDs = unique(counts[[i]]$group))
  ratio = r/mean(bardata_finetuned$Tn5Insertion)
  bardata_finetuned$downsampled_insertions = apply(bardata_finetuned['Tn5Insertion'], 1, function(x) rbinom(n = 1, size = x, p = ratio))
  Tn5_plot_finetuned = plot_insertions(bardata_finetuned$downsampled_insertions, 
                                       positions = bardata_finetuned$position, 
                                       chr = regionsBed$chr[i])
  # Tn5_plot_finetuned
  footprinting_results_finetuned = NB_footprintings(Tn5Insertion = bardata_finetuned$downsampled_insertions, 
                                                    pred_bias = bardata_finetuned$pred_bias, 
                                                    positions = bardata_finetuned$position, 
                                                    p.adjust.method = 'BY', 
                                                    nCores = 2)
  p_value_matrix_finetuned = footprinting_results_finetuned[['pval']]
  effect_size_matrix_finetuned = footprinting_results_finetuned[['effect_size']]
  
  binding_sites = binding_sites_pval(p_value_matrix_finetuned, 
                                     width_threshold = 10, 
                                     p_threshold = 0.05)
  
  if(!dir.exists(paste0("mtDNA/downsampled_comparison/", r))){
    system(paste("mkdir -p", paste0("mtDNA/downsampled_comparison/", r)))
  }
  write.table(binding_sites, 
              paste0("mtDNA/downsampled_comparison/", r, "/binding_sites_finetuned_", k, ".txt"))
  
  ### PRINT ---- 
  bardata_PRINT = get_bardata(counts = counts, 
                              regionsBed = regionsBed, 
                              bias = pred_bias_PRINT, 
                              i = i, 
                              groupIDs = unique(counts[[i]]$group))
  bardata_PRINT$downsampled_insertions = bardata_finetuned$downsampled_insertions
  Tn5_plot_PRINT = plot_insertions(bardata_PRINT$downsampled_insertions, 
                                   positions = bardata_finetuned$position, 
                                   chr = regionsBed$chr[i])
  
  footprinting_results_PRINT = NB_footprintings(Tn5Insertion = bardata_PRINT$downsampled_insertions, 
                                                pred_bias = bardata_PRINT$pred_bias, 
                                                positions = bardata_PRINT$position, 
                                                p.adjust.method = 'BH', 
                                                nCores = 2)
  p_value_matrix_PRINT = footprinting_results_PRINT[['pval']]
  effect_size_matrix_PRINT = footprinting_results_PRINT[['effect_size']]
  
 binding_sites_PRINT = binding_sites_pval(p_value_matrix_PRINT, 
                                           width_threshold = 10, 
                                           p_threshold = 0.05)
  
  if(!dir.exists(paste0("mtDNA/downsampled_comparison/", r))){
    system(paste("mkdir -p", paste0("mtDNA/downsampled_comparison/", r)))
  }
  write.table(binding_sites_PRINT, 
              paste0("mtDNA/downsampled_comparison/", r, "/binding_sites_PRINT_", k, ".txt"))
  
  
}

if (F){
  ratios = c(seq(0.1, 2, 0.1), 
             seq(2, 16, 1), 
             seq(18, 20, 2))
  
  setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/SAMN27505541"))
  
  binding_sites_finetuned = NULL
  binding_sites_PRINT = NULL
  for (r in ratios){
    for (k in 1:100){
      binding_sites = read.table(paste0("mtDNA/downsampled_comparison/", r, "/binding_sites_finetuned_", k, ".txt"), 
                                 header = TRUE)
      binding_sites_finetuned = rbind(binding_sites_finetuned, 
                                      cbind(binding_sites, 
                                            'r' = rep(r, dim(binding_sites)[1])))
      binding_sites = read.table(paste0("mtDNA/downsampled_comparison/", r, "/binding_sites_PRINT_", k, ".txt"), 
                                 header = TRUE)
      binding_sites_PRINT = rbind(binding_sites_PRINT, 
                                      cbind(binding_sites, 
                                            'r' = rep(r, dim(binding_sites)[1])))
    }
  }
  finetuned_summary = binding_sites_finetuned %>% 
    group_by(r) %>%
    summarise(count = n()/100)
  
  PRINT_summary = binding_sites_PRINT %>% 
    group_by(r) %>%
    summarise(count = n()/100)
  summary = left_join(PRINT_summary, finetuned_summary, 
                            by = 'r')
  summary[is.na(summary)] = 0
  colnames(summary)[2:3] = c('PRINT', 'TraceBIND')
  summary$PRINT[3:(dim(summary)[1] - 2)] = sapply(3:(dim(summary)[1]-2), 
                                                  function(i) mean(summary$PRINT[(i-2):(i+2)]))
  summary$TraceBIND[3:(dim(summary)[1] - 2)] = sapply(3:(dim(summary)[1]-2), 
                                                  function(i) mean(summary$TraceBIND[(i-2):(i+2)]))
  
  summary = summary %>%
    pivot_longer(cols = c('PRINT', 'TraceBIND'), 
                 names_to = 'method', 
                 values_to = 'counts')
  
  p1 = ggplot(summary, 
         aes(x = r, y = counts, color = method)) + 
    geom_point() + 
    geom_line() + 
    xlim(0, 20) + 
    labs(x = 'average reads', 
         y = 'average false positives') + 
    theme_bw()
  p1
  ggsave('false_positives_comparison_PRINT_finetuned.pdf', 
         p1, height = 3.5, width = 5)
}


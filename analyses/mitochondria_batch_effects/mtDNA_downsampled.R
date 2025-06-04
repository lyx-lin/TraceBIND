args = commandArgs(TRUE)
k = as.numeric(args[[1]])
i = 1
# reticulate::use_condaenv('base')
# projectNames = list.files('/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/CNN')
# projectName = projectNames[j]
projectNames = c('Control_6')
for (projectName in projectNames){
  setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", projectName))
  # setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/", projectName))
  source("../../../../code/utils.R")
  source("../../../../code/getBias.R")
  source("~/codes/footprints/footprint_prediction.R")
  library(tidyverse)
  # metadata = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/metadata_all.txt")
  # metadata = metadata %>%
  #   filter(library_id == projectName)
  # barcodeGroups = metadata[c('barcode', 'celltype')]
  
  counts = readRDS("mtDNA/counts2.rds")
  
  regionsBed = data.frame('chr' = 'chrM', 
                          'start' = seq(501, 16500, 1000), 
                          'end' = seq(1500, 16500, 1000))
  
  # get_bias(regionsBed, 
  #          referenceGenome = 'hg38', 
  #          path = paste0(getwd(), '/mtDNA/finetuned_comparison_1_'), 
  #          code_path = '../../../../', 
  #          model_use = paste0('Tn5_NN_model_Control_6_finetuned_mtDNA_comparison_1.h5'))
  
  fine_tuned_freezed_bias = read.table("mtDNA/finetuned_comparison_1_pred_bias.txt")
  
  # for (j in 1:9){
  #   projectName = projectNames[j]
  #   setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/", projectName))
  #   get_bias(regionsBed, 
  #            referenceGenome = 'hg38', 
  #            path = paste0(getwd(), '/mtDNA/freezed_finetuned_'), 
  #            code_path = '../../../../', 
  #            model_use = paste0('Tn5_NN_model_', projectName, '_freezed_finetuned.h5'))
  # }
  
  
  ### downsamples ----
  set.seed(k)
  ratios = seq(0.2, 3, 0.2)
  # ratios = c(seq(0.1, 5, 0.1), 
  #            seq(5.5, 40, 0.5))
  # ratios = seq(30.5, 40, 0.5)
  # if (file.exists(paste0("mtDNA/downsampled/", x, "/binding_sites_freezed_finetuned_", k, ".txt"))){
  #   next
  # }
  for (x in ratios){
    print(x)
    ### freeze finetuned
    pred_bias = fine_tuned_freezed_bias[i, ]
    # bardata$pred_bias = as.numeric(fine_tuned_freezed_bias[i, ])
    bardata = get_bardata(counts = counts, 
                          regionsBed = regionsBed, 
                          bias = pred_bias, 
                          i = i, 
                          groupIDs = unique(counts[[i]]$group))
    
    ratio = x/mean(bardata$Tn5Insertion)
    bardata$downsampled_Tn5Insertion = apply(bardata['Tn5Insertion'], 1, function(x) 
      rpois(n = 1, x*ratio))
    # bardata$Tn5Insertion = apply(bardata['Tn5Insertion'], 1, function(x) 
    #   rbinom(1, x, ratio))
    # bardata$obs_bias = bardata$Tn5Insertion 
    # footprinting_results_obsbias = NB_footprintings(Tn5Insertion = bardata$Tn5Insertion, 
    #                                                 pred_bias = bardata$obs_bias, 
    #                                                 positions = bardata$position, 
    #                                                 p.adjust.method = 'BH', 
    #                                                 nCores = 2)

    footprinting_results = NB_footprintings(Tn5Insertion = bardata$Tn5Insertion, 
                                            pred_bias = bardata$pred_bias, 
                                            positions = bardata$position, 
                                            p.adjust.method = 'BH', 
                                            nCores = 2)
    p_value_matrix = footprinting_results[['pval']]
    effect_size_matrix = footprinting_results[['effect_size']]
    
    binding_sites_freezed = binding_sites_pval(p_value_matrix, 
                                               smoothed = T, 
                                               width_threshold = 10, 
                                               p_threshold = 1)
    binding_sites_freezed['mean_coverage'] = rep(NA, dim(binding_sites_freezed)[1])

    if (dim(binding_sites_freezed)[1] > 0){
      for (index in (1:dim(binding_sites_freezed)[1])){
        positions_1 = (binding_sites_freezed$position[index] - binding_sites_freezed$width[index]/2 - 1):(binding_sites_freezed$position[index] - 1 - binding_sites_freezed$width[index]*3/2)
        positions_2 = (binding_sites_freezed$position[index] + binding_sites_freezed$width[index]/2 + 1):(binding_sites_freezed$position[index] + 1 + binding_sites_freezed$width[index]*3/2)
        mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                            mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
        
        # mean_coverages = c(mean_coverages, mean_coverage)
        # mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
        #                     mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
        binding_sites_freezed[index, 'mean_coverage'] = mean_coverage

      }
    }
    
    # binding_sites_freezed = binding_sites_pval(p_value_matrix, 
    #                                            width_threshold = 10)
    if(!dir.exists(paste0("mtDNA/downsampled_comparison_1_smoothed/", x))){
      system(paste("mkdir -p", paste0("mtDNA/downsampled_comparison_1_smoothed/", x)))
    }
    write.table(binding_sites_freezed, 
                paste0("mtDNA/downsampled_comparison_1_smoothed/", x, "/binding_sites_finetuned_", k, ".txt"))
  }
}

#### summary-----
if (F){
  # projectNames = list.files('/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/CNN')
  # projectNames = projectNames[c(10:20, 25:29)]
  projectName = 'Control_3'
  setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", 
               projectName, '/mtDNA/downsampled_freezed'))
  # setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/", 
  #             projectName, '/mtDNA/downsampled'))
  downsampled_ratios = gtools::mixedsort(list.files())
  downsampled_ratios = downsampled_ratios
  # setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA/", 
  #              projectName, '/mtDNA/downsampled_freezed/'))
  
  setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", 
               projectName, '/mtDNA/downsampled_freezed/'))
  thresholds = NULL
  for (downsampled_ratio in downsampled_ratios){
    filenames <- list.files(path = paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", 
                                          projectName, '/mtDNA/downsampled_freezed/', downsampled_ratio), 
                            pattern=".*_freezed_finetuned_.*[0-9]+", full.names=TRUE)
    filenames = gtools::mixedsort(filenames)[1:100]
    data = lapply(filenames, function(x){
      bs = read.table(x, header = TRUE) 
    }
    )
    thresholds = rbind(thresholds, 
                       data.table::rbindlist(data, fill=TRUE))
  }
  
  thresholds = as.data.frame(thresholds)
  thresholds$total_counts = thresholds$mean_coverage*thresholds$width
  # thresholds = na.omit(thresholds)
  # max(as.numeric(thresholds$mean_coverage))
  labels_total_counts = seq(0, 400, 1)
  breaks_total_counts <- c(0, 
                           sapply(1:(length(labels_total_counts) - 1), function(i) mean(c(labels_total_counts[i], labels_total_counts[i + 1]))), Inf)
  labels_total_counts <- as.character(labels_total_counts)  # Corresponding labels
  thresholds$labels_total_counts <- cut(as.numeric(thresholds$total_counts), 
                                        breaks = breaks_total_counts, 
                                        labels = labels_total_counts)
  labels_mean_counts = seq(0, 40, 0.2)
  breaks_mean_counts <- c(0, 
                          sapply(1:(length(labels_mean_counts) - 1), function(i) mean(c(labels_mean_counts[i], labels_mean_counts[i + 1]))), Inf)
  labels_mean_counts <- as.character(labels_mean_counts)  # Corresponding labels
  thresholds$labels_mean_counts <- cut(as.numeric(thresholds$mean_coverage), 
                                        breaks = breaks_mean_counts, 
                                        labels = labels_mean_counts, 
                                       right = F)
  # labels_width = seq(0, 200, 20)
  # breaks_width <- c(0, 
  #                   sapply(1:(length(labels_width) - 1), function(i) mean(c(labels_width[i], labels_width[i + 1]))), Inf)
  # labels_width <- as.character(labels_width)  # Corresponding labels
  # thresholds$labels_width <- cut(as.numeric(thresholds$width), 
  #                                breaks = breaks_width, 
  #                                labels = labels_width, 
  #                                right = F)
  
  # thresholds = na.omit(thresholds)
  logp_total_counts = thresholds %>%
    group_by(labels_total_counts) %>%
    summarise(
      count = n(), 
      # label1 = as.numeric(unique(labels)), 
      threshold = quantile(p_value, 0.9)
      # threshold = mean(p_value) + 2.32*sd(p_value)/n()
      # threshold = median(p_value) + 1.57*IQR(p_value)/sqrt(1 + n())
    )
  logp_total_counts = na.omit(logp_total_counts)
  # logp = logp[1:250, ]
  logp_total_counts$smoothed_threshold = logp_total_counts$threshold
  logp_total_counts$smoothed_threshold[3:(dim(logp_total_counts)[1] - 2)] = sapply(3:(dim(logp_total_counts)[1]-2), 
                                                         function(i) mean(logp_total_counts$threshold[(i-2):(i+2)]))
  
  logp_total_counts = logp_total_counts[1:300, ]
  
  logp_mean_counts = thresholds %>%
    group_by(labels_mean_counts) %>%
    summarise(
      count = n(), 
      # label1 = as.numeric(unique(labels)), 
      threshold = quantile(p_value, 0.9)
      # threshold = mean(p_value) + 2.32*sd(p_value)/n()
      # threshold = median(p_value) + 1.57*IQR(p_value)/sqrt(1 + n())
    )
  logp_mean_counts = na.omit(logp_mean_counts)
  # logp = logp[1:250, ]
  logp_mean_counts$smoothed_threshold = logp_mean_counts$threshold
  logp_mean_counts$smoothed_threshold[3:(dim(logp_mean_counts)[1] - 2)] = sapply(3:(dim(logp_mean_counts)[1]-2), 
                                                                                   function(i) mean(logp_mean_counts$threshold[(i-2):(i+2)]))
  logp_mean_counts = logp_mean_counts[1:150, ]
  write.table(logp_mean_counts, 
              paste0("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", projectName, "/logp_threshold_mean_counts.txt"))
  write.table(logp_total_counts, 
              paste0("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/", projectName, "/logp_threshold_total_counts.txt"))
}

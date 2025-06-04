args = commandArgs(TRUE)
k = as.numeric(args[[1]])
print(k)
i = 11
library(tidyverse)
reticulate::use_condaenv("base")
source("~/codes/footprints/footprint_prediction.R")

# projectNames = c("kidney_defnd", "kidney_multiome")
for (projectName in c("age_82_weeks/U21/kidney_multiome"# , 
                      # "age_82_weeks/U21/kidney_defnd", 
                      # "age_82_weeks/V22/kidney_defnd"# , 
                      # "age_30_weeks/L12/kidney_multiome", 
                      # "age_16_weeks/C03/kidney_multiome"
                      )){
  setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName))
  source("../../../../../code/utils.R")
  source("../../../../../code/getBias.R")

  counts = readRDS("mtDNA/counts.rds")
  
  regionsBed = data.frame('chr' = 'chrM', 
                          'start' = seq(201, 16200, 1000), 
                          'end' = seq(1200, 16200, 1000))
  # get_bias(regionsBed, 
  #          referenceGenome = 'rn7', 
  #          path = paste0(getwd(), '/mtDNA/freezed_finetuned_'), 
  #          code_path = '../../../../../', 
  #          model_use = "Tn5_NN_model_U21_multiome_freezed_finetuned.h5")
  
  fine_tuned_freezed_bias = read.table("mtDNA/freezed_finetuned_pred_bias.txt")
  
  ratios = c(seq(0.1, 5, 0.1), 
             seq(5.5, 40, 0.5))
  # ratios = seq(30.5, 30, 0.5)
  ### downsamples ----
  set.seed(k)
  for (x in ratios){
    # if (file.exists(paste0("mtDNA/downsampled/", x, "/binding_sites_freezed_finetuned_", k, ".txt"))){
    #   next
    # }
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
    bardata$downsampled_insertions = apply(bardata['Tn5Insertion'], 1, function(x) 
      rpois(n = 1, x*ratio))
    
    footprinting_results = NB_footprintings(Tn5Insertion = bardata$downsampled_insertions, 
                                            pred_bias = bardata$pred_bias, 
                                            positions = bardata$position, 
                                            p.adjust.method = 'BH', 
                                            nCores = 2)
    
    p_value_matrix = footprinting_results[['pval']]

    binding_sites_freezed = binding_sites_pval(p_value_matrix,
                                               smoothed = F, 
                                               width_threshold = 0, 
                                               p_threshold = 1)
    
    binding_sites_freezed['mean_coverage'] = rep(NA, 
                                                 dim(binding_sites_freezed)[1])
    if (dim(binding_sites_freezed)[1] > 0){
      for (index in (1:dim(binding_sites_freezed)[1])){
        positions_1 = (binding_sites_freezed$position[index] - binding_sites_freezed$width[index]/2 - 1):(binding_sites_freezed$position[index] - 1 - binding_sites_freezed$width[index]*3/2)
        positions_2 = (binding_sites_freezed$position[index] + binding_sites_freezed$width[index]/2 + 1):(binding_sites_freezed$position[index] + 1 + binding_sites_freezed$width[index]*3/2)
        # mean_coverage = min(mean(sapply(2:length(positions_1), 
        #                                 function(x) mean(bardata[bardata$position %in% positions_1[1:x], 'downsampled_insertions']))),
        #                     mean(sapply(2:length(positions_2), 
        #                                 function(x) mean(bardata[bardata$position %in% positions_2[1:x], 'downsampled_insertions']))))
        mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'downsampled_insertions']), 
                            mean(bardata[bardata$position %in% positions_2, 'downsampled_insertions']))
        
        # mean_coverages = c(mean_coverages, mean_coverage)
        # mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
        #                     mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
        binding_sites_freezed[index, 'mean_coverage'] = mean_coverage
      }
    }
    
    if(!dir.exists(paste0("mtDNA/downsampled_freezed/", x))){
      system(paste("mkdir -p", paste0("mtDNA/downsampled_freezed/", x)))
    }
    
    write.table(binding_sites_freezed, 
                paste0("mtDNA/downsampled_freezed/", x, 
                       "/binding_sites_freezed_finetuned_", k, ".txt"))
  }
}


#### summary-----
if (F){
  for (projectName in c("age_82_weeks/U21/kidney_multiome", 
                        "age_82_weeks/U21/kidney_defnd", 
                        "age_82_weeks/V22/kidney_defnd"
                        )){
    setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName, '/mtDNA/downsampled'))
    downsampled_ratios = gtools::mixedsort(list.files())
    thresholds = NULL
    for (downsampled_ratio in downsampled_ratios){
      filenames <- list.files(path = paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", 
                                            projectName, '/mtDNA/downsampled/', downsampled_ratio), 
                              pattern=".*_freezed_finetuned_.*[0-9]+", full.names=TRUE)
      filenames = gtools::mixedsort(filenames)[1:100]
      data = lapply(filenames, function(x){
        bs = read.table(x, header = TRUE) 
        bs
      })
      thresholds = rbind(thresholds, 
                         data.table::rbindlist(data, fill=TRUE))
    }
    thresholds = as.data.frame(thresholds)
    
    labels_mean_counts = seq(0.1, 40, 0.1)
    breaks_mean_counts <- c(0, 
                            sapply(1:(length(labels_mean_counts) - 1), function(i) mean(c(labels_mean_counts[i], labels_mean_counts[i + 1]))), Inf)
    labels_mean_counts <- as.character(labels_mean_counts)
    thresholds$labels_mean_counts <- cut(as.numeric(thresholds$mean_coverage), 
                                         breaks = breaks_mean_counts, 
                                         labels = labels_mean_counts)
    
    logp_mean_counts = thresholds %>%
      group_by(labels_mean_counts) %>%
      summarise(
        count = n(), 
        threshold = quantile(p_value, 0.95)
      )
    logp_mean_counts = na.omit(logp_mean_counts)

    logp_mean_counts$smoothed_threshold = logp_mean_counts$threshold
    logp_mean_counts$smoothed_threshold[3:(dim(logp_mean_counts)[1] - 2)] = sapply(3:(dim(logp_mean_counts)[1]-2), 
                                                                                   function(i) mean(logp_mean_counts$threshold[(i-2):(i+2)]))
    logp_mean_counts = logp_mean_counts[1:300, ]
    write.table(logp_mean_counts, 
                paste0("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName, "/logp_mean_counts_freezed.txt"))
    # write.table(logp_total_counts, 
    #             paste0("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName, "/logp_total_counts_freezed.txt"))
  }
}

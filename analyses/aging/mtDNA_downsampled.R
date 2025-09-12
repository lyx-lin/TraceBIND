args = commandArgs(TRUE)
k = as.numeric(args[[1]])
print(k)
i = 12

# projectNames = c("kidney_defnd", "kidney_multiome")
projectName = "L12/kidney_defnd"

setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName))
source("../../../../code/utils.R")
source("../../../../code/getBias.R")
source("~/codes/footprints/footprint_prediction.R")

# metadata = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/metadata_all.txt")
# metadata = metadata %>%
#   filter(library_id == projectName)
# barcodeGroups = metadata[c('barcode', 'celltype')]

counts = readRDS("mtDNA/counts.rds")

regionsBed = data.frame('chr' = 'chrM', 
                        'start' = seq(201, 16200, 1000), 
                        'end' = seq(1200, 16200, 1000))
# get_bias(regionsBed, 
#          referenceGenome = 'rn7', 
#          path = paste0(getwd(), '/mtDNA/freezed_finetuned_'), 
#          code_path = '../../../../', 
#          model_use = "Tn5_NN_model_L12_multiome_freezed_finetuned.h5")

fine_tuned_freezed_bias = read.table("mtDNA/freezed_finetuned_pred_bias.txt")

ratios = c(seq(0.1, 0.9, 0.2), 
           seq(1, 5, 0.1),
           seq(5.5, 15, 0.5), 
           seq(16, 30, 2))

### downsamples ----
set.seed(k)
for (x in ratios){
  if (file.exists(paste0("mtDNA/downsampled/", x, "/binding_sites_freezed_finetuned_", k, ".txt"))){
    next
  }
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
  bardata$downsampled_insertions = apply(bardata['Tn5Insertion'], 1, function(x) rbinom(n = 1, size = x, p = ratio))
  
  footprinting_results = NB_footprintings(Tn5Insertion = bardata$downsampled_insertions, 
                                          pred_bias = bardata$pred_bias, 
                                          positions = bardata$position, 
                                          p.adjust.method = 'BH', 
                                          nCores = 2)
  
  p_value_matrix = footprinting_results[['pval']]
  binding_sites_freezed = binding_sites_pval(p_value_matrix, width_threshold = 0)
  
  if(!dir.exists(paste0("mtDNA/downsampled/", x))){
    system(paste("mkdir -p", paste0("mtDNA/downsampled/", x)))
  }
  write.table(binding_sites_freezed, 
              paste0("mtDNA/downsampled/", x, "/binding_sites_freezed_finetuned_", k, ".txt"))
}

#### summary-----
if (x > 100){
  projectName = "L12/kidney_multiome"
  
  setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName, '/mtDNA/downsampled/'))
  downsampled_ratios = gtools::mixedsort(list.files())
  logp = as.data.frame(matrix(NA, 
                              ncol = length(downsampled_ratios), 
                              nrow = 1))
  colnames(logp) = downsampled_ratios
  rownames(logp) = projectName
  # projectName = projectNames[j]
  for (downsampled_ratio in downsampled_ratios){
    filenames <- list.files(path = paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", 
                                          projectName, '/mtDNA/downsampled/', downsampled_ratio), 
                            pattern=".*_freezed_finetuned_.*[0-9]+", full.names=TRUE)
    num <- unlist(lapply(filenames, function(x) max(0, (read.table(x, header = TRUE))$p_value[1], na.rm = T)))
    logp[projectName, downsampled_ratio] = mean(num)
  }
  
  
  write.table(logp, 
              paste0("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName, "/logp_threshold.txt"))
}

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Signac)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/kidney_benchmarking/NR3C1/human_kidney/Control6"))

source("~/codes/footprints/footprint_prediction.R")

counts = readRDS("../aggregated_counts.rds")

regionsBed = read.table("peaks_region.txt")
chip_positions = read.table("tf_region.txt")
chip_positions = unique(chip_positions)
logp_mean_counts = read.table("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6/logp_threshold_mean_counts.txt")

average = NULL
for (i in 1:length(counts)){
  average = c(average, sum(counts[[i]]$count)/length((regionsBed$start[i]:regionsBed$end[i])))
}
high_coverage = which(average >= 1)

## tracebind ----
process_binding_sites <- function(i) {
  file_path <- paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/results/data/NR3C1/effect_size/nb_bs_BH_freezed/", 
                      i, "_", "binding_sites_nb.txt")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  binding_sites <- read.table(file_path, header = TRUE)
  bardata <- get_bardata(counts = counts, 
                         regionsBed = regionsBed, 
                         bias = 1,
                         i = i, 
                         groupIDs = unique(counts[[i]]$group))

  binding_sites_thresholds <- data.frame(matrix(ncol = dim(binding_sites)[2], nrow = 0))
  colnames(binding_sites_thresholds) <- colnames(binding_sites)
  
  if (dim(binding_sites)[1] > 0) {
    for (index in 1:dim(binding_sites)[1]) {
      positions_1 <- (binding_sites$position[index] - binding_sites$width[index] / 2 - 1):
        (binding_sites$position[index] - 1 - binding_sites$width[index] * 3 / 2)
      positions_2 <- (binding_sites$position[index] + binding_sites$width[index] / 2 + 1):
        (binding_sites$position[index] + 1 + binding_sites$width[index] * 3 / 2)
      
      mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                          mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
      threshold <- logp_mean_counts$smoothed_threshold[(abs(mean_coverage - as.numeric(as.character(logp_mean_counts$labels_mean_counts))) == 
                                                          min(abs(mean_coverage - as.numeric(as.character(logp_mean_counts$labels_mean_counts)))))]
      
      
      if (binding_sites$p_value[index] >= mean(threshold)) {
        binding_sites_thresholds <- rbind(binding_sites_thresholds, binding_sites[index, ])
      }
    }
  }
  binding_sites = binding_sites_thresholds
  
  tf <- chip_positions[chip_positions$chr == as.character(regionsBed$chr[i]) & 
                         chip_positions$end >= regionsBed$start[i] & 
                         chip_positions$end <= regionsBed$end[i] | 
                         chip_positions$chr == as.character(regionsBed$chr[i]) & 
                         chip_positions$start >= regionsBed$start[i] & 
                         chip_positions$start <= regionsBed$end[i] | 
                         chip_positions$chr == as.character(regionsBed$chr[i]) & 
                         chip_positions$start >= regionsBed$start[i] & 
                         chip_positions$end <= regionsBed$end[i], ]
  
  tf <- unique(tf)
  tf_regions <- GenomicRanges::GRanges(seqnames = tf$chr, 
                                       ranges = IRanges::IRanges(start = tf$start, 
                                                                 end = tf$end))
  
  if (length(tf_regions) == 0) {
    return(NULL)
  }
  
  result <- list(num_chip = 0, num_footprinting_chips = 0, random_prob = c(), widths_bs = c(), num_bs = dim(binding_sites)[1])
  
  if (dim(binding_sites)[1] == 0) {
    result$num_chip <- length(tf_regions)
    return(result)
  } else {
    binding_sites_range <- data.frame("chr" = regionsBed$chr[i], 
                                      "start" = binding_sites$position_effect_size - 
                                        binding_sites$width_effect_size / 2, 
                                      "end" = binding_sites$position_effect_size + 
                                        binding_sites$width_effect_size / 2)
    
    binding_sites_regions <- GenomicRanges::GRanges(seqnames = binding_sites_range$chr, 
                                                    ranges = IRanges::IRanges(start = binding_sites_range$start, 
                                                                              end = binding_sites_range$end))
    
    width_peak <- regionsBed$end[i] - regionsBed$start[i] + 1
    result$num_chip <- length(tf_regions)
    result$num_footprinting_chips <- length(subsetByOverlaps(tf_regions, binding_sites_regions))
    result$num_bs = dim(binding_sites)[1]
    
    chip <- subsetByOverlaps(tf_regions, binding_sites_regions)
    result$random_prob <- c(result$random_prob, 
                            unlist(sapply(seq_along(chip), function(x) {
                              widths <- width(binding_sites_regions[which(countOverlaps(binding_sites_regions, chip[x]) >= 1), ])
                              widths <- mean(widths)
                              non_overlap <- max(0, start(chip)[x] - regionsBed$start[i] + 1 - widths) + 
                                max(0, regionsBed$end[i] - end(chip)[x] + 1 - widths)
                              pmin(1, 1 - non_overlap / (pmax(1, width_peak - widths)))
                            })))
    
    result$widths_bs <- c(result$widths_bs, 
                          unlist(sapply(seq_along(chip), function(x) {
                            width(GenomicRanges::reduce(binding_sites_regions[which(countOverlaps(binding_sites_regions, chip[x]) >= 1), ]))
                          })))
  }
  
  return(result)
}
results <- pbmcapply::pbmclapply(1:length(counts), 
                                 process_binding_sites, 
                                 mc.cores = 4)

num_chip <- sapply(results, function(res) if (!is.null(res)) res$num_chip else NA)
num_footprinting_chips <- sapply(results, function(res) if (!is.null(res)) res$num_footprinting_chips else NA)
random_prob <- unlist(lapply(results, function(res) if (!is.null(res)) res$random_prob else NULL))
widths_bs <- unlist(lapply(results, function(res) if (!is.null(res)) res$widths_bs else NULL))

data = readRDS("data.rds")
data[1, 1:2] = c(sum(num_footprinting_chips[high_coverage], na.rm = T)/sum(num_chip[high_coverage], na.rm = T), 
                 mean(random_prob[high_coverage], na.rm = T))
saveRDS(data, "data.rds")

### scPrinter ----
library(anndata)
tf_scprinter = read_h5ad('PRINT_TF.h5ad')
nuc_scprinter = read_h5ad('PRINT_nuc.h5ad')
tf_sites = Signac::StringToGRanges(tf_scprinter$obsm_keys(), 
                                   sep= c(":", "-"))
tf_sites = sort(tf_sites)
keys = gtools::mixedsort(tf_scprinter$obsm_keys())
regions <- GRanges(seqnames = regionsBed$chr, 
                   ranges = IRanges(start = regionsBed$start, 
                                    end = regionsBed$end))

chunkResults <- readRDS("chunkedTFBSResults/chunk_1.rds")

num_scprinter_chip = rep(0, length(counts))
num_scprinter_footprinting_chips = rep(0, length(counts))
widths_bs_scprinter = NULL
random_prob_scprinter = NULL
for (i in 1:length(keys)){
  print(i)
  regionTFBS = chunkResults[[i]]
  regionTFBS$scprinter_score = pmax(unlist(tf_scprinter$obsm[keys[i]]), 
                                    unlist(nuc_scprinter$obsm[keys[i]]))
  sites = regionTFBS$sites
  sites = sites[regionTFBS$scprinter_score >= 0.5, ]
  binding_sites = reduce(sites)
  tf = chip_positions[chip_positions$chr == as.character(regionsBed$chr[i]) & 
                        chip_positions$end >= regionsBed$start[i] & 
                        chip_positions$end <= regionsBed$end[i] | 
                        chip_positions$chr == as.character(regionsBed$chr[i]) & 
                        chip_positions$start >= regionsBed$start[i] & 
                        chip_positions$start <= regionsBed$end[i] | 
                        chip_positions$chr == as.character(regionsBed$chr[i]) & 
                        chip_positions$start >= regionsBed$start[i] & 
                        chip_positions$end <= regionsBed$end[i], ]
  
  tf = unique(tf)
  tf_regions = GenomicRanges::GRanges(seqnames = tf$chr, 
                                      ranges = IRanges::IRanges(start = tf$start, 
                                                                end = tf$end))
  
  if (length(tf_regions) == 0){
    next
  }
  if (length(binding_sites) == 0){
    num_scprinter_chip[i] = length(tf_regions)
    # num_scprinter_motifs[i] = length(motif_regions)
    num_scprinter_footprinting_chips[i] = 0
    # num_scprinter_footprinting_motifs[i] = 0
  } else {                                                                        
    # num_scprinter_motifs[i] = length(motif_regions)
    # num_scprinter_footprinting_motifs[i] = length(subsetByOverlaps(motif_regions, binding_sites))
    width_peak = regionsBed$end[i] - regionsBed$start[i] + 1
    width_chip = width(tf_regions)
    num_scprinter_chip[i] = length(tf_regions)
    num_scprinter_footprinting_chips[i] = length(subsetByOverlaps(tf_regions, binding_sites))
    chip = subsetByOverlaps(tf_regions, binding_sites)
    widths = width(reduce(binding_sites[which(countOverlaps(binding_sites, tf_regions) >= 1), ]))
    widths_bs_scprinter = c(widths_bs_scprinter,
                            unlist(sapply(seq_along(chip), function(x){
                              widths = width(reduce(binding_sites[which(countOverlaps(binding_sites, chip[x]) >= 1), ]))
                              widths[which.min(abs(widths - 16))]
                            })))
    random_prob_scprinter = c(random_prob_scprinter, 
                              unlist(sapply(seq_along(chip), function(x) {
                                widths = width(binding_sites[which(countOverlaps(binding_sites, chip[x]) >= 1), ])
                                widths = widths[which.min(abs(widths - 16))]
                                non_overlap = max(0, start(chip)[x] - regionsBed$start[i] + 1 - widths) + 
                                  max(0, regionsBed$end[i] - end(chip)[x] + 1 - widths)
                                pmin(1, 1 - non_overlap/(pmax(1, width_peak - widths)))})))
  }
}

mean(random_prob_scprinter)
sum(num_scprinter_footprinting_chips[high_coverage])/sum(num_scprinter_chip[high_coverage]) 

data = readRDS("data.rds")
data[5, 1:2] = c(sum(num_scprinter_footprinting_chips[high_coverage])/sum(num_scprinter_chip[high_coverage]), 
                 mean(random_prob_scprinter))
saveRDS(data, 'data.rds')

### HINT ----
regions = makeGRangesFromDataFrame(regionsBed)
regions = regions[high_coverage]
hint_results = read.table("HINT/footprints.bed")
colnames(hint_results)[1:3] = c('chr', 'start', 'end')
hint_results = makeGRangesFromDataFrame(hint_results)
hint_results = reduce(hint_results)
chips = makeGRangesFromDataFrame(chip_positions)
chips = unique(findOverlapPairs(chips, regions)@first)
hint_results = unique(findOverlapPairs(hint_results, chips)@first)
overlaps_HINT = findOverlaps(chips, hint_results)
length(unique(overlaps_HINT@from))/length(chips) # 0.8996798
random_prob_HINT = NULL
for (i in 1:max(overlaps_HINT@from)){
  chip = chips[i]
  if (length(chip) == 0){
    next
  }
  tf = hint_results[overlaps_HINT[overlaps_HINT@from == i]@to]
  region = findOverlapPairs(regions, chip)@first
  random_prob_HINT = c(random_prob_HINT, 
                        unlist(sapply(seq_along(chip), function(x) {
                          widths = width(tf)
                          widths = widths[which.min(abs(widths - 16))]
                          non_overlap = max(0, start(chip) - start(region) + 1 - widths) + 
                            max(0, end(region) - end(chip) + 1 - widths)
                          pmin(1, 1 - non_overlap/(pmax(1, width(region) - widths)))})))
}
mean(random_prob_HINT) # 0.2359361

### TOBIAS ----
motifs = list.dirs("TOBIAS/BINDetect_output", full.names = FALSE, recursive = FALSE)
tobias_results = data.table::rbindlist(lapply(motifs,
                                              function(x) read.table(paste0("TOBIAS/BINDetect_output/", x, '/beds/', x, '_TOBIAS_footprints_bound.bed'), sep = '\t', header = FALSE)))
colnames(tobias_results)[1:3] = c('chr', 'start', 'end')
tobias_results = makeGRangesFromDataFrame(tobias_results)
tobias_results = reduce(tobias_results)
regions = makeGRangesFromDataFrame(regionsBed)
regions = regions[high_coverage]
chips = makeGRangesFromDataFrame(chip_positions)
chips = unique(findOverlapPairs(chips, regions)@first)
tobias_results = unique(findOverlapPairs(tobias_results, chips)@first)
overlaps_TOBIAS = findOverlaps(chips, tobias_results)
length(unique(overlaps_TOBIAS@from))/length(chips) # 0.8497866
width_bs_TOBIAS = NULL
random_prob_TOBIAS = NULL
for (i in 1:max(overlaps_TOBIAS@from)){
  chip = chips[i]
  if (length(chip) == 0){
    next
  }
  tf = tobias_results[overlaps_TOBIAS[overlaps_TOBIAS@from == i]@to]
  region = findOverlapPairs(regions, chip)@first
  widths = width(tf)
  widths = widths[which.min(abs(widths - 16))]
  width_bs_TOBIAS = c(width_bs_TOBIAS, widths)
  non_overlap = max(0, start(chip) - start(region) + 1 - widths) + 
    max(0, end(region) - end(chip) + 1 - widths)
  random_prob_TOBIAS = c(random_prob_TOBIAS, 
                         pmin(1, 1 - non_overlap/(pmax(1, width(region) - widths))))
}
mean(random_prob_TOBIAS) # 0.2421622


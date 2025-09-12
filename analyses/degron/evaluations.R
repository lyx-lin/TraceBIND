library(tidyverse)
library(ggpubr)
library(GenomicRanges)
library(RColorBrewer)
library(GGally)
library(hdrcde)
library(ggdensity)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggrepel)
library(gridExtra)
source("~/codes/footprints/footprint_prediction.R")
theme_Publication <- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            # legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            # legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/degron/CTCF/"))
regionsBed = read.table("peaks_region.txt")
colnames(regionsBed) = c('chr', 'start', 'end')

counts_treatment = readRDS("treatment/counts.rds")
counts_control = readRDS("control/counts.rds")
logp_treatment = read.table("treatment/logp_threshold_freezed.txt")
logp_control = read.table("control/logp_threshold_freezed.txt")

average_treatment = NULL
for (i in 1:length(counts_treatment)){
  average_treatment = c(average_treatment, 
                        sum(counts_treatment[[i]]$count)/length(regionsBed$start[i]:regionsBed$end[i]))
}

average_control = NULL
for (i in 1:length(counts_control)){
  average_control = c(average_control, 
                      sum(counts_control[[i]]$count)/length(regionsBed$start[i]:regionsBed$end[i]))
}

motifs_regions = read.table('motifs_regions.txt')
### negative and positive ----
chip_control = read.table("chip_region.txt")
# chip_control = chip_control[chip_control$V7 >= fprcutoff, ]
chip_range = makeGRangesFromDataFrame(chip_control)
overlaps = unique(findOverlaps(motifs_ranges, chip_range)@from)
motifs_regions = motifs_regions[overlaps, ]
motifs_ranges = makeGRangesFromDataFrame(motifs_regions)
motifs_ranges = reduce(motifs_ranges)

chip_treatment_regions = read.table("treatment/ENCFF881BGO.bed")
colnames(chip_treatment_regions)[1:3] = c('chr', 'start', 'end')
chip_treatment_ranges = makeGRangesFromDataFrame(chip_treatment_regions)
overlaps_treatment_control = unique(findOverlaps(chip_treatment_ranges, motifs_ranges))
length(overlaps_treatment_control)

motifs_regions$motif_treatment = 0
motifs_regions$motif_treatment[overlaps_treatment_control@to] = overlaps_treatment_control@from

### TraceBIND ------
motifs_regions$peak_index = sapply(1:length(motifs_ranges), function(i) which(regionsBed$chr == as.character(motifs_regions$seqnames[i]) & 
                regionsBed$end <= motifs_regions$start[i] & 
                regionsBed$end >= motifs_regions$end[i] | 
                regionsBed$chr == as.character(motifs_regions$seqnames[i]) & 
                regionsBed$start <= motifs_regions$start[i] & 
                regionsBed$start >= motifs_regions$end[i] | 
                regionsBed$chr == as.character(motifs_regions$seqnames[i]) & 
                regionsBed$start <= motifs_regions$start[i] & 
                regionsBed$end >= motifs_regions$end[i], )[1])
motifs_regions = unique(motifs_regions)
motifs_regions = motifs_regions[!is.na(motifs_regions$peak_index), ]
motifs_ranges = makeGRangesFromDataFrame(motifs_regions)

process_motifs = function(i){
  index = motifs_regions$peak_index[i]
  effect_size_matrix_control = read.table(paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/CTCF/control_freezed/", 
                                                 index, "_effect_size_matrix.txt"))
  p_value_matrix_control = read.table(paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/CTCF/control_freezed/", 
                                             index, "_p_value_matrix.txt"))
  colnames(p_value_matrix_control) = gsub("X", "", colnames(p_value_matrix_control))
  log_p_value_control = -log10(p_value_matrix_control)
  log_p_value_control = as.matrix(log_p_value_control)
  
  rownames(log_p_value_control) = rownames(p_value_matrix_control)
  colnames(log_p_value_control) = colnames(p_value_matrix_control)
  effect_size_matrix_control = apply(effect_size_matrix_control, 2, caTools::runmax, 5)
  effect_size_matrix_control = apply(effect_size_matrix_control, 2, pracma::conv, 5)
  effect_size_matrix_control = effect_size_matrix_control/(2*5)
  rownames(effect_size_matrix_control) = rownames(p_value_matrix_control)
  colnames(effect_size_matrix_control) = colnames(p_value_matrix_control)
  ranges = max(regionsBed$start[index], motifs_regions$start[i]):min(regionsBed$end[index], motifs_regions$end[i])
  
  bardata <- get_bardata(counts = counts_control, 
                         regionsBed = regionsBed, 
                         bias = rep(1, length(regionsBed$start[index]:regionsBed$end[index])),
                         i = index, 
                         groupIDs = unique(counts_control[[index]]$group))
  
  positions_1 <- min(regionsBed$start[index], motifs_regions$end[i] + 1):min(regionsBed$end[index], motifs_regions$end[i] + 30)
  positions_2 <- max(regionsBed$start[index], motifs_regions$start[i] - 30):max(regionsBed$start[index], motifs_regions$start[i] - 1)
  
  mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                      mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
  
  threshold <- logp_control$smoothed_threshold[which.min(abs(mean_coverage - as.numeric(rownames(logp_control))))]
  
  control_log10_pval = max(log_p_value_control[as.character(ranges), as.character(2*(5:25))])
  control_effect_size = ifelse(control_log10_pval >= max(1, threshold),  
                                 max(effect_size_matrix_control[as.character(ranges), as.character(2*(5:25))]), 
                                 min(effect_size_matrix_control[as.character(ranges), as.character(2*(5:25))]))
  
  effect_size_matrix_treatment = read.table(paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/CTCF/treatment_freezed/", 
                                                   index, "_effect_size_matrix.txt"))
  p_value_matrix_treatment = read.table(paste0("~/nzhanglab/project/linyx/footprints/results/data/degron/CTCF/treatment_freezed/", 
                                               index, "_p_value_matrix.txt"))
  colnames(p_value_matrix_treatment) = gsub("X", "", colnames(p_value_matrix_treatment))
  log_p_value_treatment = -log10(p_value_matrix_treatment)
  log_p_value_treatment = as.matrix(log_p_value_treatment)
  
  rownames(log_p_value_treatment) = rownames(p_value_matrix_treatment)
  colnames(log_p_value_treatment) = colnames(p_value_matrix_treatment)
  effect_size_matrix_treatment = apply(effect_size_matrix_treatment, 2, caTools::runmax, 5)
  effect_size_matrix_treatment = apply(effect_size_matrix_treatment, 2, pracma::conv, 5)
  effect_size_matrix_treatment = effect_size_matrix_treatment/(2*5)
  rownames(effect_size_matrix_treatment) = rownames(p_value_matrix_treatment)
  colnames(effect_size_matrix_treatment) = colnames(p_value_matrix_treatment)
  ranges = max(regionsBed$start[index], motifs_regions$start[i]):min(regionsBed$end[index], motifs_regions$end[i])
  
  bardata <- get_bardata(counts = counts_treatment, 
                         regionsBed = regionsBed, 
                         bias = rep(1, length(regionsBed$start[index]:regionsBed$end[index])),
                         i = index, 
                         groupIDs = unique(counts_treatment[[index]]$group))
  
  positions_1 <- min(regionsBed$start[index], motifs_regions$end[i] + 1):min(regionsBed$end[index], motifs_regions$end[i] + 30)
  positions_2 <- max(regionsBed$start[index], motifs_regions$start[i] - 30):max(regionsBed$start[index], motifs_regions$start[i] - 1)
  
  mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'Tn5Insertion']), 
                      mean(bardata[bardata$position %in% positions_2, 'Tn5Insertion']))
  
  threshold <- logp_treatment$smoothed_threshold[which.min(abs(mean_coverage - as.numeric(rownames(logp_treatment))))]
  
  treatment_log10_pval = max(log_p_value_treatment[as.character(ranges), as.character(2*(5:25))])
  treatment_effect_size = ifelse(treatment_log10_pval >= max(1, threshold),  
                                 max(effect_size_matrix_treatment[as.character(ranges), as.character(2*(5:25))]), 
                                 min(effect_size_matrix_treatment[as.character(ranges), as.character(2*(5:25))]))
  result <- list(control_log10_pval = control_log10_pval, 
                 control_effect_size = control_effect_size, 
                 treatment_log10_pval = treatment_log10_pval, 
                 treatment_effect_size = treatment_effect_size
  )
  return(result)
}

results = pbmcapply::pbmclapply(1:dim(motifs_regions)[1], 
                                process_motifs, 
                                mc.cores = 2)
motifs_regions$control_log10_pval = sapply(results, function(res) res$control_log10_pval)
motifs_regions$control_effect_size = sapply(results, function(res) res$control_effect_size)
motifs_regions$treatment_log10_pval = sapply(results, function(res) res$treatment_log10_pval)
motifs_regions$treatment_effect_size = sapply(results, function(res) res$treatment_effect_size)


plot(PRROC::pr.curve(motifs_regions$treatment_log10_pval1[motifs_regions$motif_treatment == 1], 
                     motifs_regions$treatment_log10_pval1[motifs_regions$motif_treatment == 0], 
                     curve = T))

plot(PRROC::pr.curve(as.numeric(motifs_regions$treatment_effect_size1[motifs_regions$motif_treatment == 1]), 
                     as.numeric(motifs_regions$treatment_effect_size1[motifs_regions$motif_treatment == 0]), 
                     curve = T))

write.table(motifs_regions, "motifs_regions_score.txt")

### scPRINT ---- 
library(anndata)
tf_control_scprinter = read_h5ad('control/PRINT_TF.h5ad')
tf_treatment_scprinter = read_h5ad('treatment/PRINT_TF.h5ad')

tf_control_sites = Signac::StringToGRanges(tf_control_scprinter$obsm_keys(), 
                                            sep= c(":", "-"))
tf_control_sites = sort(tf_control_sites)
keys_control = gtools::mixedsort(tf_control_scprinter$obsm_keys())
tf_treatment_sites = Signac::StringToGRanges(tf_treatment_scprinter$obsm_keys(), 
                                            sep= c(":", "-"))
tf_treatment_sites = sort(tf_treatment_sites)
keys_treatment = gtools::mixedsort(tf_treatment_scprinter$obsm_keys())

regions <- GRanges(seqnames = regionsBed$chr, 
                   ranges = IRanges(start = regionsBed$start, 
                                    end = regionsBed$end))
chunkResults_control = readRDS("control/chunkedTFBSResults/chunk_1.rds")
chunkResults_treatment = readRDS("treatment/chunkedTFBSResults/chunk_1.rds")

process_motifs_scPrinter = function(i){
  index = motifs_regions$peak_index[i]
  regionTFBS_control = chunkResults_control[[index]]
  regionTFBS_control$scprinter_score = unlist(tf_control_scprinter$obsm[keys_control[index]])
  regionTFBS_treatment = chunkResults_treatment[[index]]
  regionTFBS_treatment$scprinter_score = unlist(tf_treatment_scprinter$obsm[keys_treatment[index]])
  
  sites = regionTFBS_control$sites
  ranges = motifs_ranges[i]
  
  j = findOverlaps(sites, ranges)@from
  if (length(j) > 0){
    score_control = max(regionTFBS_control$scprinter_score[j])
    score_treatment = max(regionTFBS_treatment$scprinter_score[j])
  } else {
    score_control = 0
    score_treatment = 0
  }
  
  return(list(scprinter_control_score = score_control, 
              scprinter_treatment_score = score_treatment))
}
scprinter_results = pbmcapply::pbmclapply(1:length(motifs_ranges), 
                                          process_motifs_scPrinter, 
                                          mc.cores = 2)
motifs_regions$scprinter_control = sapply(scprinter_results, function(res) res$scprinter_control_score)
motifs_regions$scprinter_treatment = sapply(scprinter_results, function(res) res$scprinter_treatment_score)
plot(PRROC::pr.curve(motifs_regions$scprinter_treatment[motifs_regions$motif_treatment > 0], 
                     motifs_regions$scprinter_treatment[motifs_regions$motif_treatment == 0], 
                     curve = T))

### PRINT ----
chunkResults_control = readRDS("control/chunkedTFBSResults/chunk_1.rds")
chunkResults_treatment = readRDS("treatment/chunkedTFBSResults/chunk_1.rds")

process_motifs_PRINT = function(i){
  index = motifs_regions$peak_index[i]
  regionTFBS_control = chunkResults_control[[index]]
  regionTFBS_treatment = chunkResults_treatment[[index]]
  
  sites = regionTFBS_control$sites
  ranges = motifs_ranges[i]
  
  j = findOverlaps(sites, ranges)@from
  if (length(j) > 0){
    score_control = max(regionTFBS_control$TFBSScores[j])
    score_treatment = max(regionTFBS_treatment$TFBSScores[j])
  } else {
    score_control = 0
    score_treatment = 0
  }
  
  return(list(PRINT_control_score = score_control, 
              PRINT_treatment_score = score_treatment))
}
motifs_regions$PRINT_control_score = PRINT_control_score
motifs_regions$PRINT_treatment_score = PRINT_treatment_score

write.table(motifs_regions, 
            'motif_region_score.txt')

PRINT_results = pbmcapply::pbmclapply(1:length(motifs_ranges), 
                                      process_motifs_PRINT, 
                                      mc.cores = 2)

### HINT ----
hint_control_results = read.table("control/HINT/footprints.bed")
hint_treatment_results = read.table("treatment/HINT/footprints.bed")

colnames(hint_control_results)[1:3] = c('chr', 'start', 'end')
colnames(hint_treatment_results)[1:3] = c('chr', 'start', 'end')

hint_control_results = makeGRangesFromDataFrame(hint_control_results, 
                                                keep.extra.columns = T)
hint_treatment_results = makeGRangesFromDataFrame(hint_treatment_results, 
                                                  keep.extra.columns = T)

process_motifs_HINT = function(i){
  ranges = motifs_ranges[i]
  j = findOverlaps(hint_control_results, ranges)@from
  if (length(j) > 0){
    score_control = max(hint_control_results[j]$V5)
  } else {
    score_control = 0
  }
  j = findOverlaps(hint_treatment_results, ranges)@from
  if (length(j) > 0){
    score_treatment = max(hint_treatment_results[j]$V5)
  } else {
    score_treatment = 0
  }
  return(list(score_treatment = score_treatment, 
              score_control = score_control))
}
hint_results = pbmcapply::pbmclapply(1:length(motifs_ranges), 
                                     process_motifs_HINT, 
                                     mc.cores = 2)
motifs_regions$control_hint = sapply(hint_results, function(res) res$score_control)
motifs_regions$treatment_hint = sapply(hint_results, function(res) res$score_treatment)

### Tobias ----
tobias_control_results = read.table("control/tobias/CTCF_MA0139.2/beds/CTCF_MA0139.2_all.bed")
tobias_treatment_results = read.table("treatment/tobias/CTCF_MA0139.2/beds/CTCF_MA0139.2_all.bed")

colnames(tobias_control_results)[1:3] = c('chr', 'start', 'end')
colnames(tobias_treatment_results)[1:3] = c('chr', 'start', 'end')

tobias_control_results = makeGRangesFromDataFrame(tobias_control_results, 
                                                keep.extra.columns = T)
tobias_treatment_results = makeGRangesFromDataFrame(tobias_treatment_results, 
                                                  keep.extra.columns = T)

process_motifs_tobias = function(i){
  ranges = motifs_ranges[i]
  j = findOverlaps(tobias_control_results, ranges)@from
  if (length(j) > 0){
    score_control = max(tobias_control_results[j]$V10)
  } else {
    score_control = 0
  }
  j = findOverlaps(tobias_treatment_results, ranges)@from
  if (length(j) > 0){
    score_treatment = max(tobias_treatment_results[j]$V10)
  } else {
    score_treatment = 0
  }
  return(list(score_treatment = score_treatment, 
              score_control = score_control))
}
tobias_results = pbmcapply::pbmclapply(1:length(motifs_ranges), 
                                     process_motifs_tobias, 
                                     mc.cores = 2)
motifs_regions$control_tobias = sapply(tobias_results, function(res) res$score_control)
motifs_regions$treatment_tobias = sapply(tobias_results, function(res) res$score_treatment)
plot(PRROC::pr.curve(motifs_regions$treatment_tobias[motifs_regions$motif_treatment == 1], 
                     motifs_regions$treatment_tobias[motifs_regions$motif_treatment == 0], 
                     curve = T))

### roc and pr -----
motifs_regions = read.table('motifs_regions_score.txt')
motifs_ranges = makeGRangesFromDataFrame(motifs_regions)

motifs_regions_negative = motifs_regions[motifs_regions$motif_treatment == 0, ]
motifs_regions_positive = NULL

for (i in unique(motifs_regions$motif_treatment[motifs_regions$motif_treatment != 0])){
  motifs_regions_positive = rbind(motifs_regions_positive, 
                                  cbind(motifs_regions[which(motifs_regions$motif_treatment == i)[1], 1:3], 
                                        t(apply(motifs_regions[motifs_regions$motif_treatment == i, 
                                                               4:18], 2, max))))
}
motifs_regions_positive = data.frame(motifs_regions_positive)
motifs_regions = rbind(motifs_regions_negative, motifs_regions_positive)
dim(motifs_regions)

motifs_regions[, 2:18] <- lapply(motifs_regions[, 2:18], as.numeric)

motifs_regions = motifs_regions[motifs_regions$peak_index %in% which(average_treatment > 0.15 & average_control > 0.15), ]
sum(motifs_regions$motif_treatment == 0)
sum(motifs_regions$motif_treatment != 0)

log10_treatment_pr = PRROC::pr.curve(motifs_regions$treatment_log10_pval[motifs_regions$motif_treatment != 0], 
                                     motifs_regions$treatment_log10_pval[motifs_regions$motif_treatment == 0], 
                                     curve = T)

enrichment_treatment_pr = PRROC::pr.curve(motifs_regions$treatment_effect_size[motifs_regions$motif_treatment != 0], 
                                          motifs_regions$treatment_effect_size[motifs_regions$motif_treatment == 0], 
                                          curve = T)
enrichment_treatment_pr$auc.integral
PRINT_treatment_pr = PRROC::pr.curve(motifs_regions$PRINT_treatment_score[motifs_regions$motif_treatment != 0], 
                                     motifs_regions$PRINT_treatment_score[motifs_regions$motif_treatment == 0], 
                                     curve = T)
scprinter_treatment_pr = PRROC::pr.curve(motifs_regions$scprinter_treatment[motifs_regions$motif_treatment != 0], 
                                         motifs_regions$scprinter_treatment[motifs_regions$motif_treatment == 0], 
                                         curve = T)

hint_treatment_pr = PRROC::pr.curve(motifs_regions$treatment_hint[motifs_regions$motif_treatment != 0], 
                                    motifs_regions$treatment_hint[motifs_regions$motif_treatment == 0], 
                                    curve = T)
tobias_treatment_pr = PRROC::pr.curve(motifs_regions$treatment_tobias[motifs_regions$motif_treatment != 0], 
                                      motifs_regions$treatment_tobias[motifs_regions$motif_treatment == 0], 
                                      curve = T)

# Convert to data frame for ggplot
log10_treatment_pr_df <- data.frame(Recall = log10_treatment_pr$curve[, 1], 
                                    Precision = log10_treatment_pr$curve[, 2], 
                                    method = paste0("NBinomial -log(pval)\n(AUC = ", 
                                                    round(log10_treatment_pr$auc.integral, 2), ")"))

enrichment_treatment_pr_df <- data.frame(Recall = enrichment_treatment_pr$curve[, 1], 
                                         Precision = enrichment_treatment_pr$curve[, 2], 
                                         method = paste0('TraceBind, AUC = ',
                                                         round(enrichment_treatment_pr$auc.integral, 2)))

PRINT_treatment_pr_df <- data.frame(Recall = PRINT_treatment_pr$curve[, 1], 
                                    Precision = PRINT_treatment_pr$curve[, 2], 
                                    method = paste0('PRINT, AUC = ', 
                                                    round(PRINT_treatment_pr$auc.integral, 2)))
scprinter_treatment_pr_df <- data.frame(Recall = scprinter_treatment_pr$curve[, 1], 
                                        Precision = scprinter_treatment_pr$curve[, 2], 
                                        method = paste0('seq2PRINT, AUC = ', 
                                                        round(scprinter_treatment_pr$auc.integral, 2)))

hint_treatment_pr_df <- data.frame(Recall = hint_treatment_pr$curve[, 1], 
                                   Precision = hint_treatment_pr$curve[, 2], 
                                   method = paste0('HINT, AUC = ', 
                                                   round(hint_treatment_pr$auc.integral, 2)))


tobias_treatment_pr_df <- data.frame(Recall = tobias_treatment_pr$curve[, 1], 
                                     Precision = tobias_treatment_pr$curve[, 2], 
                                     method = paste0('TOBIAS, AUC = ', 
                                                     round(tobias_treatment_pr$auc.integral, 2)))


pr_df = rbind(# log10_treatment_pr_df, 
  enrichment_treatment_pr_df, 
  scprinter_treatment_pr_df, 
  PRINT_treatment_pr_df, 
  hint_treatment_pr_df, 
  tobias_treatment_pr_df)
pr_df$method <- factor(pr_df$method, 
                       levels = c(paste0('TraceBind, AUC = ',
                                round(enrichment_treatment_pr$auc.integral, 2)), 
                         paste0('seq2PRINT, AUC = ', 
                                round(scprinter_treatment_pr$auc.integral, 2)), 
                         paste0('PRINT, AUC = ', 
                                round(PRINT_treatment_pr$auc.integral, 2)), 
                         paste0('HINT, AUC = ', 
                                round(hint_treatment_pr$auc.integral, 2)), 
                         paste0('TOBIAS, AUC = ', 
                                round(tobias_treatment_pr$auc.integral, 2))))# Plot Precision-Recall Curve using ggplot2
pr = ggplot(pr_df, aes(x = Recall, y = Precision)) +
  geom_line(aes(color = method), linewidth = 2) +
  labs(x = "Recall",
       y = "Precision") +
  theme_bw() + 
  ylim(0.45, 1) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8, 'Pastel1')[c(1:3, 4:5)]) + 
  theme_Publication(base_size = 12) + 
  theme(legend.position = c(0.7, 0.2), 
        legend.title=element_blank())

pr 

ggsave("pr_auc.pdf", pr, 
       height = 5, width = 5.5)
ggsave("pr_auc.png", pr, 
       height = 5, width = 5.5)

### roc  ----
roc_enrich_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                  motifs_regions$treatment_effect_size)
roc_logp_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                motifs_regions$treatment_log10_pval)
roc_PRINT_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                 motifs_regions$PRINT_treatment_score)
roc_scprint_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                   motifs_regions$scprinter_treatment)
roc_hint_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                   motifs_regions$treatment_hint)
roc_tobias_treatment <- pROC::roc(motifs_regions$motif_treatment, 
                                   motifs_regions$treatment_tobias)

# Convert ROC curves to data frames
roc_enrich_treatment_df <- data.frame(FPR = 1 - (roc_enrich_treatment$specificities), 
                                      TPR = roc_enrich_treatment$sensitivities, 
                                      method = paste("NBinomial_enrich_score\n(AUC =", 
                                                     round(roc_enrich_treatment$auc, 2), ")"))

roc_logp_treatment_df <- data.frame(FPR = 1- (roc_logp_treatment$specificities), 
                                    TPR = roc_logp_treatment$sensitivities, 
                                    method = paste("NBinomial_logp_pval\n(AUC =", 
                                                   round(roc_logp_treatment$auc, 2), ")"))

roc_PRINT_treatment_df <- data.frame(FPR = 1 - (roc_PRINT_treatment$specificities), 
                                     TPR = roc_PRINT_treatment$sensitivities, 
                                     method = paste("PRINT score\n(AUC =", 
                                                    round(roc_PRINT_treatment$auc, 2), ")"))
roc_scprint_treatment_df <- data.frame(FPR = 1 - (roc_scprint_treatment$specificities), 
                                     TPR = roc_scprint_treatment$sensitivities, 
                                     method = paste("seq2PRINT score\n(AUC =", 
                                                    round(roc_scprint_treatment$auc, 2), ")"))
roc_hint_treatment_df <- data.frame(FPR = 1 - (roc_hint_treatment$specificities), 
                                     TPR = roc_hint_treatment$sensitivities, 
                                     method = paste("HINT score\n(AUC =", 
                                                    round(roc_hint_treatment$auc, 2), ")"))
roc_tobias_treatment_df <- data.frame(FPR = 1 - (roc_tobias_treatment$specificities), 
                                    TPR = roc_tobias_treatment$sensitivities, 
                                    method = paste("TOBIAS score\n(AUC =", 
                                                   round(roc_tobias_treatment$auc, 2), ")"))


# Combine both ROC curves into one data frame
roc_treatment <- rbind(roc_enrich_treatment_df, 
                       # roc_logp_treatment_df, 
                       # roc_PRINT_treatment_df, 
                       roc_scprint_treatment_df, 
                       roc_hint_treatment_df, 
                       roc_tobias_treatment_df
                       )

roc_treatment$method <- factor(roc_treatment$method, 
                               levels = c(paste("NBinomial_enrich_score\n(AUC =", 
                                       round(roc_enrich_treatment$auc, 2), ")"), 
                                 paste("seq2PRINT score\n(AUC =", 
                                       round(roc_scprint_treatment$auc, 2), ")"), 
                                 # paste0('PRINT score\n(AUC = ', 
                                 #        round(PRINT_treatment_pr$auc.integral, 2), ")"), 
                                 paste("HINT score\n(AUC =", 
                                       round(roc_hint_treatment$auc, 2), ")"), 
                                 paste("TOBIAS score\n(AUC =", 
                                       round(roc_tobias_treatment$auc, 2), ")")))# Plot Precision-Recall Curve using ggplot2

roc = ggplot(roc_treatment, 
             aes(x = FPR, y = TPR)) +
  geom_line(aes(color = method), linewidth = 2) +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, 'Pastel1')[c(1:2, 5, 4)]) + 
  theme_Publication() 
roc

ggsave("roc_auc.png", roc, 
       height = 5, width = 8.5)
ggsave("roc_auc.pdf", roc, 
       height = 5, width = 8.5)

### before and after ====
motifs_undepleted = motifs_regions[motifs_regions$motif_treatment != 0, ]
motifs_depleted = motifs_regions[motifs_regions$motif_treatment == 0, ]

motifs_beforeafter <- motifs_regions[motifs_regions$motif_treatment == 0, ] %>%
  mutate(id = row_number()) %>%  # Add row identifier for tracking
  pivot_longer(
    cols = c(control_effect_size, 
             treatment_effect_size, 
             control_log10_pval, 
             treatment_log10_pval, 
             # PRINT_control_score, 
             # PRINT_treatment_score, 
             scprinter_control, 
             scprinter_treatment, 
             control_hint, 
             treatment_hint, 
             control_tobias,
             treatment_tobias),
    names_to = "method_group",
    values_to = "score"
  ) %>%
  mutate(
    group = case_when(
      grepl("control", method_group) ~ "control",
      grepl("treatment", method_group) ~ "treatment"
    ),
    method = case_when(
      grepl("effect_size", method_group) ~ "TraceBind",
      grepl("scprinter", method_group) ~ "seq2PRINT", 
      grepl("hint", method_group) ~ "HINT", 
      grepl("tobias", method_group) ~ "TOBIAS"
    )
  ) %>%
  dplyr::select(group, score, method)  %>%
  group_by(method) %>%
  mutate(score_scaled = ifelse(method == c('TraceBind'), 
                               score, score / quantile(score, 1, na.rm = TRUE))) %>%
  
  ungroup()
n = sum(motifs_regions$motif_treatment == 0)
mu <- n*(n+1)/4
sigma <- sqrt(n*(n+1)*(2*n + 1)/24)
wilcox_tracebind = ks.test(motifs_regions[motifs_regions$motif_treatment == 0, 'control_effect_size'], 
                               motifs_regions[motifs_regions$motif_treatment == 0, 'treatment_effect_size'])$statistic
wilcox_stat_tracebind = (wilcox_tracebind-mu)/sigma

wilcox_hint = ks.test(motifs_regions[motifs_regions$motif_treatment == 0, 'control_hint'], 
                               motifs_regions[motifs_regions$motif_treatment == 0, 'treatment_hint'])$statistic
wilcox_stat_hint = (wilcox_hint-mu)/sigma

wilcox_tobias = ks.test(motifs_regions[motifs_regions$motif_treatment == 0, 'control_tobias'], 
                               motifs_regions[motifs_regions$motif_treatment == 0, 'treatment_tobias'])$statistic
wilcox_stat_tobias = (wilcox_tobias-mu)/sigma

wilcox_scprinter = ks.test(motifs_regions[motifs_regions$motif_treatment == 0, 'scprinter_control'], 
                            motifs_regions[motifs_regions$motif_treatment == 0, 'scprinter_treatment'])$statistic
wilcox_stat_scprinter = (wilcox_scprinter-mu)/sigma

wilcox_stat = data.frame('method' = c('HINT', 'seq2PRINT', 'TOBIAS', 'TraceBind'), 
                         'wilcox_stat' = c(wilcox_stat_hint, 
                                           wilcox_stat_scprinter, 
                                           wilcox_stat_tobias, 
                                           wilcox_stat_tracebind))

before_after = ggplot(data = motifs_beforeafter) + 
  geom_boxplot(aes(x = method, y = score_scaled, fill = group), 
               outlier.shape = 16, outlier.size = 1, alpha = 0.8) + 
  geom_text(data = wilcox_stat, 
            aes(x = method, y = 1.1, 
                label = paste0("U stat = ", round(wilcox_stat, 0))), 
            inherit.aes = FALSE, size = 4) + 
  scale_x_discrete(limits = c('TraceBind', 
                              "seq2PRINT", 
                              "HINT", "TOBIAS")) + 
  scale_fill_manual(values = c("#FFCB42", "#99DBF5")) + 
  scale_y_continuous(limits = c(-0.5, 1.1), 
                     breaks = seq(-0.5, 1, 0.25)) + 
  labs(x = "Method", y = "Scaled Score", fill = "Group") + 
  theme_Publication(base_size = 18)
before_after
ggsave("before_after_logfc.png", before_after, 
       height = 5, width = 6)
ggsave("before_after_logfc.pdf", before_after, 
       height = 5, width = 6)

projectName = 'CTCF'
print(projectName)
setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/degron/", projectName))

library(MASS)
library(tidyverse)
library(ggplot2)
library(GGally)
library(patchwork)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

peaks_control = read.table("control/outs/filtered_peak_bc_matrix/peaks.bed")
peaks_treatment = read.table("treatment/outs/filtered_peak_bc_matrix/peaks.bed")
colnames(peaks_control) = c('chr', 'start', 'end')
colnames(peaks_treatment) = c('chr', 'start', 'end')

peaks_control = makeGRangesFromDataFrame(peaks_control)
peaks_treatment = makeGRangesFromDataFrame(peaks_treatment)
overlaps = findOverlapPairs(peaks_control, peaks_treatment)
peaks = pbmcapply::pbmclapply(1:length(overlaps), function(i) 
  reduce(c(overlaps@first[i], overlaps@second[i])), 
  mc.cores = 2
)
peaks = as(peaks, "GRangesList")
peaks = unlist(peaks)
peaks = peaks[seqnames(peaks) %in% paste0('chr', 1:22)]

write.table(peaks, "peaks_region.txt", sep = '\t', )

peaks = read.table("peaks_region.txt")
peaks = makeGRangesFromDataFrame(peaks)

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection="CORE",tax_group='vertebrates',all_versions=FALSE)
)
tf = pwm['MA0139.1']

chip_control = read.table("control/ENCFF463FGL.bed")
colnames(chip_control)[1:3] = c('chr', 'start', 'end')
chip_control_granges = makeGRangesFromDataFrame(chip_control, 
                                                keep.extra.columns = T)
chip_control_granges = chip_control_granges[unique(findOverlaps(chip_control_granges, peaks)@from)]

motif_control_pos = matchMotifs(tf, 
                                chip_control_granges, 
                                genome = "hg38", 
                                out = "positions")

motif_control_pos = makeGRangesFromDataFrame(motif_control_pos)
chip_overlaps = findOverlaps(motif_control_pos, 
                             chip_control_granges)
chip_control_granges = chip_control_granges[unique(chip_overlaps@to)]
motif_granges = motif_control_pos[unique(chip_overlaps@from)]

peaks_overlap = findOverlaps(motif_control_pos, peaks)
peaks = peaks[unique(peaks_overlap@to)]
motif_pos = unique(findOverlapPairs(motif_control_pos, 
                                    peaks)@first)
regionsBed = data.frame("chr" = seqnames(peaks), 
                        "start" = start(peaks), 
                        "end" = end(peaks))
write.table(regionsBed, "peaks_region.txt", sep = '\t', )
write.table(motif_control_pos, "motif_region.txt", sep = '\t')
write.table(as.data.frame(chip_control_granges), "chip_region.txt", sep = '\t')

metadata = read.table("control/outs/filtered_peak_bc_matrix/barcodes.tsv")
metadata$group = 'HCT116'
colnames(metadata) = c('barcodeID', 'group')
rownames(metadata) = NULL

get_bias(regionsBed, 
         referenceGenome = 'hg38', 
         path = paste0(getwd(), '/freezed_finetuned_'), 
         code_path = '../../../../', 
         model_use = 'Tn5_NN_model_CTCF_control_freezed_finetuned.h5'
         )
frags_path = paste0("control/outs/fragments.tsv.gz")
counts = get_count(frags_path, 
                   regionsBed = regionsBed, 
                   barcodeGroups = metadata, 
                   chunkSize = 6000)
saveRDS(counts, 'control/counts.rds')

metadata = read.table("treatment/outs/filtered_peak_bc_matrix/barcodes.tsv")
metadata$group = 'HCT116'
colnames(metadata) = c('barcodeID', 'group')
rownames(metadata) = NULL

get_bias(regionsBed, 
         referenceGenome = 'hg38', 
         path = paste0(getwd(), '/freezed_finetuned_'), 
         code_path = '../../../../', 
         model_use = 'Tn5_NN_model_CTCF_treatment_freezed_finetuned.h5'
)
frags_path = paste0("treatment/outs/fragments.tsv.gz")
counts = get_count(frags_path, 
                   regionsBed = regionsBed, 
                   barcodeGroups = metadata, 
                   chunkSize = 6000)
saveRDS(counts, 'treatment/counts.rds')
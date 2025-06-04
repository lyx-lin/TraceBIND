args <- commandArgs(TRUE)
i = as.numeric(args[[1]])
set.seed(1)

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
model_names = c("U21_multiome", 
                "U21_defnd", 
                "V22_defnd")

projectName = projectNames[i]
model_name = model_names[i]

print(c(projectName, model_name))
setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/", projectName))

source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("../../peaks_region.txt")
colnames(regionsBed) = c('chr', 'start', 'end')

sample_name = strsplit(model_name, '_')[[1]][1]
tech_name = strsplit(model_name, '_')[[1]][2]
tech_name = ifelse(tech_name == 'multiome', 'Multiome', 'DEFND')
barcodeGroups = read.csv(paste0("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-", sample_name, "-KYC-LN-01-", tech_name, "-ATAC/barcodes.tsv"), 
                         header = FALSE)
colnames(barcodeGroups) = 'barcodeID'
barcodeGroups$group = '1'
barcodeGroups = barcodeGroups[sample(1:dim(barcodeGroups)[1], 3500), ]
frags_path = paste0("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-", sample_name, "-KYC-LN-01-", tech_name, "-ATAC/fragments.tsv.gz")
counts = get_count(frags_path, 
                   regionsBed = regionsBed, 
                   barcodeGroups = barcodeGroups, 
                   chunkSize = 6000)
saveRDS(counts, 'counts_subset.rds')

fasta_file <- readDNAStringSet("~/nzhanglab/data/ParkerWilson_DEFND/refernece/fasta/genome.fa", format = "fasta")
genome_list <- setNames(as.list(fasta_file), names(fasta_file))
names(fasta_file) = sapply(names(fasta_file), function(x) strsplit(x, ' ')[[1]][1])

regions = GenomicRanges::GRanges(seqnames = regionsBed$chr, 
                                 ranges = IRanges::IRanges(start = regionsBed$start, 
                                                           end = regionsBed$end))
contextRadius = 50
contextLen = 2*contextRadius + 1
index = which(end(regions) + contextLen > width(fasta_file[seqnames(regions)]))
regions = regions[-index]
regionsBed = regionsBed[-index, ]
write.table(regionsBed, 
            "../../peaks_region.bed")
regionSeqs = NULL
for (regionInd in 1:length(regions)){
  range = regions[regionInd]
  regionSeq = subseq(fasta_file[[seqnames(range)]], 
                     start = start(range) - contextRadius, 
                     width = width(range) + contextLen - 1)
  regionSeqs = c(regionSeqs, 
                 as.character(regionSeq)
  )
}

regionSeqs = as.character(regionSeqs)
write.table(regionSeqs, 
            paste("freezed_finetuned_regionSeqs.txt", sep = ""), quote = F,
            col.names = F, row.names = F)
get_bias(regionsBed, 
         referenceGenome = 'hg38', 
         path = paste0(getwd(), '/freezed_finetuned_'), 
         code_path = '../../../../../', 
         model_use = paste0('Tn5_NN_model_', model_name, '_freezed_finetuned.h5')
         )


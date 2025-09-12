args <- commandArgs(TRUE)
i = as.numeric(args[[1]])
reticulate::use_condaenv("base")

library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reticulate)
library(tensorflow)
library(Biostrings)
library(GenomicRanges)

set.seed(1112)
projectNames = c('age_82_weeks/U21/kidney_multiome', 
                 'age_56_weeks/Q17/kidney_multiome', 
                 'age_30_weeks/L12/kidney_multiome', 
                 'age_16_weeks/C03/kidney_multiome')
projectName = projectNames[i]
print(c(projectName))

sample_name = strsplit(strsplit(projectName, '_')[[1]][3], '/')[[1]][2]

setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/", projectName))

source("~/codes/footprints/footprint_prediction.R")
regionsBed = read.table("../../../all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')

fasta_file <- readDNAStringSet("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/refernece/fasta/genome.fa", format = "fasta")
genome_list <- setNames(as.list(fasta_file), names(fasta_file))
names(fasta_file) = sapply(names(fasta_file), function(x) strsplit(x, ' ')[[1]][1])

regions = GenomicRanges::GRanges(seqnames = regionsBed$chr, 
                                 ranges = IRanges::IRanges(start = regionsBed$start, 
                                                           end = regionsBed$end))
contextRadius = 50
contextLen = 2*contextRadius + 1
index = which(end(regions) + contextLen > width(fasta_file[seqnames(regions)]))
if (length(index) > 0){
  regions = regions[-index]
  regionsBed = regionsBed[-index, ]
}

write.table(regionsBed, 
            "../../../all_peaks_region.bed")
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

metadata = read.table("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/metadata.txt")
metadata = metadata[metadata$library_type == 'Multiome',]
metadata = metadata[metadata$library_id == paste0("RAGE24-", sample_name, "-KYC-LN-01-Multiome-RNA"), ]
metadata = metadata['celltype']
metadata$GEX_bc = rownames(metadata)
metadata$GEX_bc = gsub(substr(metadata$GEX_bc[1], nchar(metadata$GEX_bc[1]) + 1 - 2, nchar(metadata$GEX_bc[1])), '', metadata$GEX_bc)

cells_conversion = read.csv("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/barcode_conversion_10Xmultiome.csv")
barcodeGroups = left_join(metadata, cells_conversion, by = 'GEX_bc')
barcodeGroups$GEX_bc = NULL
barcodeGroups = barcodeGroups[c(2, 1)]
barcodeGroups$ATAC_bc = paste0(barcodeGroups$ATAC_bc, '-1')
barcodeGroups = barcodeGroups[barcodeGroups$celltype %in% c("PT", "PT-MT"), ]

frags_path = paste0("~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-", sample_name, "-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz")
counts = get_count(frags_path, 
                   regionsBed = regionsBed, 
                   barcodeGroups = barcodeGroups, 
                   chunkSize = 5000)
saveRDS(counts, 'counts.rds')

barcodeGroups = barcodeGroups[sample(1:dim(barcodeGroups)[1], 2000), ]
counts_subset = get_count(frags_path, 
                   regionsBed = regionsBed, 
                   barcodeGroups = barcodeGroups, 
                   chunkSize = 5000)
saveRDS(counts_subset, 'counts_subset.rds')

get_bias(regionsBed, 
         path = paste0(getwd(), '/freezed_finetuned_'), 
         code_path = '../../../../../../', 
         model_use = paste0('Tn5_NN_model_', sample_name, '_multiome_freezed_finetuned.h5')
)
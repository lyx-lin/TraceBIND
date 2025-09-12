reticulate::use_condaenv("base")

library(MASS)
library(tidyverse)
library(ggplot2)
library(GGally)
library(patchwork)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6"))
source("/home/mnt/weka/nzh/team/woodsqu2/codes/footprints/github/models.R")
source("/home/mnt/weka/nzh/team/woodsqu2/codes/footprints/github/utils.R")
source("/home/mnt/weka/nzh/team/woodsqu2/codes/footprints/github/visualizations.R")

barcodes = read.csv("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/data/ParkerWilson/kidney_multiome_atlas/kidney_chasm/Control_6/bcanno.csv")
barcodes = barcodes[c('barcode', 'celltype')]

regionsBed = data.frame('chr' = 'chrM', 
                        'start' = seq(501, 16500, 1000), 
                        'end' = seq(1500, 16500, 1000))

### obtain counts for each cell type and each base pair 
counts = get_count("mtDNA/fragments_collapse_within.bed", 
                   regionsBed, 
                   barcodes)
saveRDS(counts, 'counts.rds')
counts = readRDS('counts.rds')

regionsBed = data.frame('chr' = 'chrM', 
                        'start' = 501:16500, 
                        'end' = 501:16500)
regions = makeGRangesFromDataFrame(regionsBed)
ranges = IRanges::IRanges(start = regionsBed$start, 
                          end = regionsBed$end)

### obtain DNA sequence 
contextRadius = 50
contextLen = 2*contextRadius + 1
regionSeqs = pbmcapply::pbmclapply(
  1:length(regions),
  function(regionInd){
    range = regions[regionInd]
    regionSeq = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                                   as.character(seqnames(range)), 
                                   start = start(range) - contextRadius, 
                                   width = width(range) + contextLen - 1,
                                   as.character = T)
    regionSeq
  },
  mc.cores = 2
)

context = do.call(rbind.data.frame, regionSeqs)
index_n = which(apply(context, 1, function(x) grepl('N', x, fixed = TRUE)))

insertions = NULL
for (i in 1:16){
  insertions = c(insertions, get_insertions(countData = counts, 
                                            regionInd = i, 
                                            groupIDs = unique(counts[[i]]$group), 
                                            width = 1000))
}

### get relative Tn5 preference
obsbias = data.frame('context' = context, 
                     'insertions' = as.numeric(insertions))
obsbias$positions = 501:16500
colnames(obsbias) = c('context', 'insertions', 'positions')

obsbias$insertions = as.numeric(obsbias$insertions)
obsbias$positions = as.numeric(obsbias$positions)
obsbias['obs_bias'] = apply(obsbias, 1, function(x){
  as.numeric(x[2])/as.numeric(mean(obsbias[as.numeric(obsbias$positions) >= as.numeric(x[3]) - 50 & 
                                             as.numeric(obsbias$positions) <= as.numeric(x[3]) + 50, "insertions"]))
}
)

rownames(obsbias) = NULL
obsbias$BACInd = rep(1:16, each = 1000)
write.table(obsbias, 'obsBias.tsv', sep = "\t")
obsBias_finetuned = obsbias[-index_n, ]
write.table(obsBias_finetuned, 'obsBias_finetuned.tsv', sep = "\t")

### Finetuned the model 
finetuned_model(code_path = "/home/mnt/weka/nzh/team/woodsqu2/codes/footprints/github/", 
                obsbias_path = "/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6/mtDNA/obsBias_finetuned.tsv", 
                PRINT_model_path = "/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/shared/Tn5_NN_model.h5", 
                finetuned_model_save_path = "/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/data/shared", 
                finetuned_model_name = "Tn5_NN_model_Control_6_freezed_finetuned"
                )


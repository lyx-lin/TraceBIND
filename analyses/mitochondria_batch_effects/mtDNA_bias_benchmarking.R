library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(tidyverse)
library(ggplot2)
library(GGally)
library(hdrcde)
library(ggdensity)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggpubr)
library(ggrepel)
library(gridExtra)
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
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

source("~/codes/footprints/footprint_prediction.R")
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
index = index
peaks = makeGRangesFromDataFrame(data.frame('chr' = 'chrM', 
                                            'start' = seq(501, 16500, 1000), 
                                            'end' = seq(1500, 16500, 1000)))
regionsBed = data.frame('chr' = 'chrM', 
                        'start' = seq(501, 16500, 1000), 
                        'end' = seq(1500, 16500, 1000))
### finetuned and PRINT ----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
files = list.files()[c(1:3, 5, 13:17)]
for (file in files){
  setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
  setwd(paste0(file))
  get_bias(regionsBed, 
           referenceGenome = 'hg38', 
           path = paste0(getwd(), '/mtDNA/finetuned_'), 
           code_path = '../../../../', 
           model_use = paste0('Tn5_NN_model_', file, '_finetuned.h5'))
  # get_bias(regionsBed, 
  #          referenceGenome = 'hg38', 
  #          path = paste0(getwd(), '/mtDNA/PRINT_'), 
  #          code_path = '../../../../', 
  #          model_use = 'Tn5_NN_model.h5')
}

setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
files = list.files()[c(8:18, 20:24, 28:33)]
for (file in files){
  setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
  setwd(paste0(file))
  get_bias(regionsBed, 
           referenceGenome = 'hg38', 
           path = paste0(getwd(), '/mtDNA/finetuned_'), 
           code_path = '../../../../', 
           model_use = paste0('Tn5_NN_model_', file, '_finetuned.h5'))
  # get_bias(regionsBed, 
  #          referenceGenome = 'hg38', 
  #          path = paste0(getwd(), '/mtDNA/PRINT_'), 
  #          code_path = '../../../../')
}

onehotEncode <- function(seq){
  onehotMapping <- list("A" = 1, "C" = 2, "G" = 3, "T" = 4, "N" = 0)
  len <- nchar(seq)
  onehotInd <- sapply(
    substring(seq, 1:len , 1:len),
    function(x){onehotMapping[[x]]})
  onehot <- matrix(0, nrow = 4, ncol = len)
  undeterminedBases <- onehotInd == 0
  onehot[cbind(unname(onehotInd), 1:len)[!undeterminedBases,]] <- 1
  onehot
}
PWMScoring <- function(seq, PWM){
  onehotSeq <- onehotEncode(seq)
  sum(PWM * onehotSeq)
}

### PWM multiome -----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
files = list.files()[c(1:3, 5, 13:17)]

radius = c(10, 15, 25)
for (PWMRadius in radius){
  PWMsBias = matrix(NA, nrow = 1000, ncol = length(files))
  rownames(PWMsBias) = 1:1000
  colnames(PWMsBias) = files
  # cor = data.frame(matrix(NA, nrow = length(files), ncol = 1))
  # rownames(cor) = files
  # colnames(cor) = 'PWM'
  for (j in c(1:9)){
    counts = readRDS(paste0(files[j], "/mtDNA/counts.rds"))
    obsbias = read.table(paste0(files[j], "/mtDNA/obsBias.tsv"))
    obsbias = obsbias$obs_bias[obsbias$BACInd == 11]
    tilePFM <- pbmcapply::pbmclapply(
      c(1:10, 12:16),
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        ATACTrack <- get_insertions(counts, 
                                    i, 
                                    unique(counts[[i]]$group), 
                                    width(peaks[i]))
        
        # One-hot encoding of each sequence context that appeared in the BAC
        # This is used to quantify the background sequence frequencies
        backgroundOnehot <- lapply(
          1:length(contextSeq), 
          function(index){onehotEncode(contextSeq[index])})
        backgroundPFM <- Reduce("+", backgroundOnehot)
        
        # Since each context is cut with different frequencies, we calculate
        # the foreground frequencies by weighing the background with cutting density 
        foregroundOnehot <- lapply(
          1:length(contextSeq), 
          function(index){backgroundOnehot[[index]] * ATACTrack[index]})
        foregroundPFM <- Reduce("+", foregroundOnehot)
        
        rownames(foregroundPFM) <- c("A", "C", "G", "T")
        rownames(backgroundPFM) <- c("A", "C", "G", "T")
        
        list(foregroundPFM, backgroundPFM)
      }, 
      mc.cores = 2
    )
    fgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[1]]}))
    bgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[2]]}))
    background <- rowSums(bgPFM)
    background <- background / sum(background)
    
    # Get PPM 
    PPM <- t(t(fgPFM) / colSums(fgPFM))
    adjustedPPM <- PPM / background
    PWM <- log2(adjustedPPM)
    
    PWMsBias[, j] = unlist(pbmcapply::pbmclapply(
      11,
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        PWMBias  <-  2 ^ sapply(contextSeq, function(seq){PWMScoring(seq, PWM)})
        PWMBias
      }, 
      mc.cores = 2
    ))
  }
  saveRDS(PWMsBias, 
          paste0(PWMRadius, '_PWMsBias.rds'))
}

PWMsBias = readRDS("25_PWMsBias.rds")
p1 = GGally::ggpairs(PWMsBias)
ggsave("corrlation_PWMsBias_multiome_25.png", 
       p1, height = 20, width = 20)

### PWM snATAC -----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
files = list.files()[c(13:18, 21:25, 34:38, 42:47)]

radius = c(10, 15, 25)
for (PWMRadius in radius){
  PWMsBias = matrix(NA, nrow = 1000, ncol = length(files))
  rownames(PWMsBias) = 1:1000
  colnames(PWMsBias) = files
  
  for (j in 1:length(files)){
    setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
    counts = readRDS(paste0(files[j], "/mtDNA/counts.rds"))
    tilePFM <- pbmcapply::pbmclapply(
      c(1:10, 12:16),
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        ATACTrack <- get_insertions(counts, 
                                    i, 
                                    unique(counts[[i]]$group), 
                                    width(peaks[i]))
        
        # One-hot encoding of each sequence context that appeared in the BAC
        # This is used to quantify the background sequence frequencies
        backgroundOnehot <- lapply(
          1:length(contextSeq), 
          function(index){onehotEncode(contextSeq[index])})
        backgroundPFM <- Reduce("+", backgroundOnehot)
        
        # Since each context is cut with different frequencies, we calculate
        # the foreground frequencies by weighing the background with cutting density 
        foregroundOnehot <- lapply(
          1:length(contextSeq), 
          function(index){backgroundOnehot[[index]] * ATACTrack[index]})
        foregroundPFM <- Reduce("+", foregroundOnehot)
        
        rownames(foregroundPFM) <- c("A", "C", "G", "T")
        rownames(backgroundPFM) <- c("A", "C", "G", "T")
        
        list(foregroundPFM, backgroundPFM)
      }, 
      mc.cores = 2
    )
    fgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[1]]}))
    bgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[2]]}))
    background <- rowSums(bgPFM)
    background <- background / sum(background)
    
    # Get PPM 
    PPM <- t(t(fgPFM) / colSums(fgPFM))
    adjustedPPM <- PPM / background
    PWM <- log2(adjustedPPM)
    
    PWMsBias[, j] = unlist(pbmcapply::pbmclapply(
      11,
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        PWMBias  <-  2 ^ sapply(contextSeq, function(seq){PWMScoring(seq, PWM)})
        PWMBias
      }, 
      mc.cores = 2
    ))
  }
  saveRDS(PWMsBias, 
          paste0(PWMRadius, '_PWMsBias.rds'))
  
}

PWMsBias = readRDS('25_PWMsBias.rds')
p1 = GGally::ggpairs(PWMsBias)

ggsave("corrlation_PWMsBias_snATAC_25.png", 
       p1, height = 20, width = 20)

### kmer ----
# multiome 
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
files = list.files()[c(1:3, 5, 13:17)]
kmerBiasMatrix = matrix(NA, nrow = 1000, ncol = length(files))
rownames(kmerBiasMatrix) = 1:1000
colnames(kmerBiasMatrix) = files

radius = c(3, 5, 7)
for (k in radius){
  for(j in 1:length(files)){
    counts = readRDS(paste0(files[j], "/mtDNA/counts.rds"))
    obsbias = read.table(paste0(files[j], "/mtDNA/obsBias.tsv"))
    bases <- c("A", "C", "G", "T")
    kmers <- as.matrix(eval(parse(text = paste("expand.grid(", paste(rep("bases", k), collapse = ", "), ")"))))
    kmers <- sapply(1:dim(kmers)[1], function(x){paste(kmers[x,], collapse = "")})
    
    # Go through each genomic tile and record insertion numbers for kmers
    tileKmer <- pbmcapply::pbmclapply(
      c(1:10, 12:16),
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = k)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        ATACTrack <- get_insertions(counts, 
                                    i, 
                                    unique(counts[[i]]$group), 
                                    width(peaks[i]))
        kmerFG <- rep(0, length(kmers))
        names(kmerFG) <- kmers
        kmerBG <- rep(0, length(kmers))
        names(kmerBG) <- kmers
        kmerFG[contextSeq] <- ATACTrack
        kmerBG[contextSeq] <- 1
        list(kmerFG, kmerBG)
      }
    )
    kmerFG <- colSums(bind_rows(sapply(tileKmer, function(x){x[[1]]})), na.rm = T)
    kmerBG <- colSums(bind_rows(sapply(tileKmer, function(x){x[[2]]})), na.rm = T)
    kmerBias <- kmerFG / kmerBG
    
    kmerBiasMatrix[, j] = unlist(pbmcapply::pbmclapply(
      11,
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = k)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        kmerBias <-  unname(kmerBias[contextSeq])
        kmerBias
      }, 
      mc.cores = 2
    ))
  }
  
  saveRDS(kmerBiasMatrix, paste0(k, "mer_bias.rds"))
}

# snATAC
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
files = list.files()[c(13:18, 21:25, 34:38, 42:47)]
kmerBiasMatrix = matrix(NA, nrow = 1000, ncol = length(files))
rownames(kmerBiasMatrix) = 1:1000
colnames(kmerBiasMatrix) = files

radius = c(3, 5, 7)

for (k in radius){
  for(j in 1:length(files)){
    counts = readRDS(paste0(files[j], "/mtDNA/counts.rds"))
    obsbias = read.table(paste0(files[j], "/mtDNA/obsBias.tsv"))
    bases <- c("A", "C", "G", "T")
    kmers <- as.matrix(eval(parse(text = paste("expand.grid(", paste(rep("bases", k), collapse = ", "), ")"))))
    kmers <- sapply(1:dim(kmers)[1], function(x){paste(kmers[x,], collapse = "")})
    
    # Go through each genomic tile and record insertion numbers for kmers
    tileKmer <- pbmcapply::pbmclapply(
      c(1:10, 12:16),
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = k)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        ATACTrack <- get_insertions(counts, 
                                    i, 
                                    unique(counts[[i]]$group), 
                                    width(peaks[i]))
        kmerFG <- rep(0, length(kmers))
        names(kmerFG) <- kmers
        kmerBG <- rep(0, length(kmers))
        names(kmerBG) <- kmers
        kmerFG[contextSeq] <- ATACTrack
        kmerBG[contextSeq] <- 1
        list(kmerFG, kmerBG)
      }
    )
    kmerFG <- colSums(bind_rows(sapply(tileKmer, function(x){x[[1]]})), na.rm = T)
    kmerBG <- colSums(bind_rows(sapply(tileKmer, function(x){x[[2]]})), na.rm = T)
    kmerBias <- kmerFG / kmerBG
    
    kmerBiasMatrix[, j] = unlist(pbmcapply::pbmclapply(
      11,
      function(i){
        tileRange <- peaks[i]
        tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
        contextRanges <- IRanges::resize(tilePositions, fix = "center", width = k)
        contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
        kmerBias <-  unname(kmerBias[contextSeq])
        kmerBias
      }, 
      mc.cores = 2
    ))
  }
  
  saveRDS(kmerBiasMatrix, paste0(k, "mer_bias.rds"))
}

cor(kmerBiasMatrix, 
    use = 'complete.obs')

kmerBiasMatrix = readRDS('5mer_bias.rds')
ggplot(data = kmerBiasMatrix, 
       aes(x = Control_4, 
           y = SAMN27505542)) +
  ylab("sample 1") + 
  xlab("sample 1") + 
  geom_hdr(aes(fill = after_stat(probs)), 
           color = "black", alpha = 0.8, 
           probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", 
           label.x = -0.5, 
           digits = 3, 
           label.y = 4000) +
  xlim(c(-0.5, 7500)) +
  # ylim(c(-0.5, 15)) +
  # ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()
ggsave("../7_mer.pdf", height = 5, width = 5
)

### plot ----
snATAC_kmer = readRDS("../snATAC_mtDNA/5mer_bias.rds")
snATAC_kmer = snATAC_kmer[, c("Control_1", "Control_2", "Control_3",
                            "Control_4", "Control_5", "Control_6",
                            "DN_1", "DN_2", "DN_3", "DN_4", "DN_5", 
                            "CKD_1", "CKD_2", "CKD_3", "CKD_4", "CKD_5")]
p1 = ggplot() +
  geom_point(aes(x = snATAC_kmer[, 'Control_1'], 
                 y = snATAC_kmer[, 'Control_4'])) + 
  xlab("Control_1") + 
  ylab("Control_4") + 
  theme_classic()
p1
ggsave("../snATAC_mtDNA/correlation_5-mer_control_1_4.png", 
       p1, height = 4, width = 4)
p2 = GGally::ggpairs(as.data.frame(snATAC_kmer))
ggsave("../snATAC_mtDNA/correlation_5_mer.png", 
       p2, height = 12, width = 12)

## benchmarking ----
cor_data = as.data.frame(matrix(ncol = 6, nrow = 0))
colnames(cor_data) = c('mito-finetuned', 'PRINT', 
                       'PWM (r=10)', 'PWM (r=25)', 
                       'kmer (k=3)', 'kmer (k=5)')
### multiome
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
files = list.files()[c(1:3, 5, 13:17)]

setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
PWMsBias_25 = readRDS('25_PWMsBias.rds')
PWMsBias_10 = readRDS('10_PWMsBias.rds')
kmer_3 = readRDS('3mer_bias.rds')
kmer_5 = readRDS('5mer_bias.rds')
kmer_7 = readRDS('7mer_bias.rds')

for (file in files){
  setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/multiome_kidney_mtDNA")
  setwd(paste0(file))
  obs_bias = read.table('mtDNA/obsBias.tsv')
  obs_bias = obs_bias[obs_bias$BACInd == 11, ]
  print_bias = read.table('mtDNA/PRINT_pred_bias.txt')
  print_bias = as.numeric(print_bias[11, ])
  finetuned_freezed_bias = read.table('mtDNA/finetuned_pred_bias.txt')
  finetuned_freezed_bias = as.numeric(finetuned_freezed_bias[11, ])
  cor_data[file, ] = c(cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + finetuned_freezed_bias)), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + print_bias)), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(PWMsBias_10[, file])), 
                           use = 'complete.obs'), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(PWMsBias_25[, file])), 
                           use = 'complete.obs'), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(kmer_3[, file])), 
                           use = 'complete.obs'), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(kmer_5[, file])), 
                           use = 'complete.obs')
                       )
}

### snATAC
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
files = list.files()[c(8:18, 20:24, 28:33)]

PWMsBias_50 = readRDS("50_PWMsBias_benchmarking.rds")
PWMsBias_25 = readRDS("25_PWMsBias_benchmarking.rds")
PWMsBias_10 = readRDS("10_PWMsBias_benchmarking.rds")
kmer_3 = readRDS('3mer_bias_benchmarking.rds')
kmer_5 = readRDS('5mer_bias_benchmarking.rds')
kmer_7 = readRDS('7mer_bias_benchmarking.rds')

for (file in files){
  setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA")
  setwd(paste0(file))
  obs_bias = read.table('mtDNA/obsBias.tsv')
  obs_bias = obs_bias[obs_bias$BACInd == 11, ]
  print_bias = read.table('mtDNA/PRINT_pred_bias.txt')
  print_bias = as.numeric(print_bias[11, ])
  finetuned_freezed_bias = read.table('mtDNA/finetuned_pred_bias.txt')
  finetuned_freezed_bias = as.numeric(finetuned_freezed_bias[11, ])
  
  cor_data[file, ] = c(cor(log(0.01 + obs_bias$obs_bias), 
                          log(0.01 + finetuned_freezed_bias)), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + print_bias)), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(PWMsBias_25[, file]))), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(PWMsBias_10[, file]))), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(kmer_3[, file])), 
                           use = 'complete.obs'), 
                       cor(log(0.01 + obs_bias$obs_bias), 
                           log(0.01 + as.numeric(kmer_5[, file])), 
                           use = 'complete.obs')
                      )
}

cor_data_visualiztion = cor_data %>%
  tibble::rownames_to_column("sample") %>%
  pivot_longer(c('mito-finetuned', 'PRINT', 
                 'PWM (r=25)', 'PWM (r=10)', 
                 'kmer (k=3)', 'kmer (k=5)'
                 ), 
               names_to = "method", 
               values_to = "correlation")
cor_data_visualiztion$method = factor(cor_data_visualiztion$method, 
                         levels = c('mito-finetuned', 'PRINT', 
                                    'PWM (r=25)', 'PWM (r=10)', 
                                    'kmer (k=5)', 'kmer (k=3)'
                                    ))
p1 = ggplot(data = cor_data_visualiztion, 
       aes(x = method, 
           y = correlation, 
           fill = method)) + 
  geom_boxplot(alpha=0.8) + 
  # geom_jitter(width = 0.2) + 
  # geom_text(aes(label=cor, y = cor), 
  #           position=position_dodge(width=0.9), vjust=0, size=3) +
  ylab('correlation with observed Tn5 bias') + 
  xlab("Methods") + 
  ylim(c(0.35, 1)) +
  scale_fill_manual(values = c("#FBB4AE", 
                               '#B3CDE3', 
                               '#DDEB9D', 
                               '#DDEB9D', 
                               '#FDEB9D', 
                               '#FDEB9D'
                               ), 
                    guide="none") +
  theme_Publication()
p1
ggsave("../../adjusted_log_correlation_box_all.pdf", 
       height = 4, width = 7)

## obsbias, PRINT, finetuned -----
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/batch_effects/snATAC_mtDNA/Control_6")

obsbias = read.table("mtDNA/obsBias_finetuned.tsv")
PRINT_bias = read.table("mtDNA/PRINT_pred_bias.txt")
finetuned_bias = read.table("mtDNA/finetuned_pred_bias.txt")

data = data.frame('obs_bias' = obsbias$obs_bias[obsbias$BACInd == 11], 
                  'PRINT_bias' = as.numeric(PRINT_bias[11, ]), 
                  'finetuned_bias' = as.numeric(finetuned_bias[11, ]))
p1 = ggplot(data = data, 
       aes(x = log(0.01 + PRINT_bias), 
           y = log(0.01 + obs_bias))) + 
  geom_point(size = 0.8) + 
  # geom_hdr(aes(fill = after_stat(probs)), 
  #          color = "black", alpha = 0.8, 
  #          probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  ylab('observed bias') +
  xlab('PRINT bias') + 
  # scale_fill_manual(values = brewer.pal(5, "Reds")) +
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  theme_Publication()

p2 = ggplot(data = data, 
            aes(x = log(0.01 + finetuned_bias), 
                y = log(0.01 + obs_bias))) + 
  geom_point(size = 0.8) + 
  # geom_hdr(aes(fill = after_stat(probs)), 
  #          color = "black", alpha = 0.8, 
  #          probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  ylab('observed bias') +
  xlab('mito-finetuned bias') + 
  # scale_fill_manual(values = brewer.pal(5, "Pastel")) +
  theme_Publication()
p1/p2

ggsave('../comparison_PRINT_finetuned_bias.pdf', 
       p1/p2, 
       height = 6, width = 3.5)

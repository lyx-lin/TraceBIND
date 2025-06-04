library(chromVAR)
library(Signac)
library(Seurat)
library(tidyverse)
library(GenomicRanges)
library(motifmatchr)
library(JASPAR2020)
library(Matrix)
library(TFBSTools)
setwd("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")
regionsBed = read.table("all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')
peaks = makeGRangesFromDataFrame(regionsBed)

metadata = read.table("~/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/metadata.txt")
cells_conversion = read.csv("~/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/barcode_conversion_10Xmultiome.csv")

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection="CORE",
              tax_group='vertebrates',
              all_versions=FALSE)
)

genome_file <- Rsamtools::FaFile("~/team/woodsqu2/nzhanglab/data/ParkerWilson_DEFND/reference/fasta/genome.fa")

motif_pos = matchMotifs(pwm, 
                        subject = peaks, 
                        genome = genome_file)

motif_matrix = as.matrix(motif_pos@assays@data$motifMatches)
library(Matrix)
motif_matrix <- Matrix::Matrix(motif_matrix, sparse = TRUE)
motif_matrix <- as(motif_matrix, "dMatrix")
colnames(motif_matrix) = names(pwm)
rownames(motif_matrix) = paste0(regionsBed$chr, '-', 
                                regionsBed$start, '-', 
                                regionsBed$end)

saveRDS(motif_matrix, 'motif_matrix.rds')
motif_matrix = readRDS("motif_matrix.rds")

## binding sites ----
binding_sites_C03 = readRDS('age_16_weeks/C03/thresholded_binding_sites.rds')
binding_sites_L12 = readRDS('age_30_weeks/L12/thresholded_binding_sites.rds')
binding_sites_Q17 = readRDS('age_56_weeks/Q17/thresholded_binding_sites.rds')
binding_sites_U21 = readRDS('age_82_weeks/U21/thresholded_binding_sites.rds')

binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill = T)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill = T)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill = T)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill = T)

binding_sites_C03$start = binding_sites_C03$position_effect_size - binding_sites_C03$width_effect_size/2
binding_sites_C03$end = binding_sites_C03$position_effect_size + binding_sites_C03$width_effect_size/2
binding_sites_L12$start = binding_sites_L12$position_effect_size - binding_sites_L12$width_effect_size/2
binding_sites_L12$end = binding_sites_L12$position_effect_size + binding_sites_L12$width_effect_size/2
binding_sites_Q17$start = binding_sites_Q17$position_effect_size - binding_sites_Q17$width_effect_size/2
binding_sites_Q17$end = binding_sites_Q17$position_effect_size + binding_sites_Q17$width_effect_size/2
binding_sites_U21$start = binding_sites_U21$position_effect_size - binding_sites_U21$width_effect_size/2
binding_sites_U21$end = binding_sites_U21$position_effect_size + binding_sites_U21$width_effect_size/2

binding_sites_C03_ranges = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size <= 120, ], 
                                                    keep.extra.columns = T)
binding_sites_L12_ranges = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size <= 120, ], 
                                                    keep.extra.columns = T)
binding_sites_Q17_ranges = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size <= 120, ], 
                                                    keep.extra.columns = T)
binding_sites_U21_ranges = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size <= 120, ], 
                                                    keep.extra.columns = T)

### seuratobject  ----
seuratobject = readRDS("seuratobject_whole_peaks.rds")
# seuratobject = seuratobject %>%
#   subset(celltype %in% c("PT", "PT-MT"))
motif_pos = matchMotifs(pwm, 
                        subject = peaks, 
                        genome = genome_file, 
                        out = 'position')
### 16 weeks =====
motif_fp_matrix = pbmcapply::pbmclapply(seq_along(motif_pos), 
                                        function(i){
                                          motif_peaks = rep(0, length(peaks))
                                          motifs = motif_pos[[i]]
                                          overlaps = findOverlaps(motifs, binding_sites_C03_ranges)
                                          motifs = motifs[overlaps@from]
                                          motif_peaks[unique(findOverlaps(peaks, motifs)@from)] = 1
                                          return(motif_peaks)
                                        }, mc.cores = 2)
motif_fp_matrix <- do.call(cbind, motif_fp_matrix)
motif_fp_matrix = Matrix::Matrix(motif_fp_matrix, sparse = TRUE)
colnames(motif_fp_matrix) = names(pwm)
rownames(motif_fp_matrix) = paste0(regionsBed$chr, '-', 
                                   regionsBed$start, '-', 
                                   regionsBed$end)
weeks16 = seuratobject %>%
  subset(dataset == 'weeks_16')
weeks16 <- RunChromVAR(
  object = weeks16,
  assay = 'ATAC',
  motif.matrix = motif_fp_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks16[['chromvar']], "age_16_weeks/C03/chromvar_footprints_vertebrates_0.05_120.rds")

weeks16 <- RunChromVAR(
  object = weeks16,
  assay = 'ATAC',
  motif.matrix = motif_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks16[['chromvar']], "age_16_weeks/C03/chromvar_vertebrates.rds")

### 30 weeks =====
motif_fp_matrix = pbmcapply::pbmclapply(seq_along(motif_pos), 
                                        function(i){
                                          motif_peaks = rep(0, length(peaks))
                                          motifs = motif_pos[[i]]
                                          overlaps = findOverlaps(motifs, binding_sites_L12_ranges)
                                          motifs = motifs[overlaps@from]
                                          motif_peaks[unique(findOverlaps(peaks, motifs)@from)] = 1
                                          return(motif_peaks)
                                        }, mc.cores = 2)
motif_fp_matrix <- do.call(cbind, motif_fp_matrix)
motif_fp_matrix = Matrix::Matrix(motif_fp_matrix, sparse = TRUE)
colnames(motif_fp_matrix) = names(pwm)
rownames(motif_fp_matrix) = paste0(regionsBed$chr, '-', 
                                   regionsBed$start, '-', 
                                   regionsBed$end)
weeks30 = seuratobject %>%
  subset(dataset == 'weeks_30')
weeks30 <- RunChromVAR(
  object = weeks30,
  assay = 'ATAC',
  motif.matrix = motif_fp_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks30[['chromvar']], "age_30_weeks/L12/chromvar_footprints_vertebrates_0.05_120.rds")

weeks30 <- RunChromVAR(
  object = weeks30,
  assay = 'ATAC',
  motif.matrix = motif_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks30[['chromvar']], "age_30_weeks/L12/chromvar_vertebrates.rds")

### 56 weeks =====
motif_fp_matrix = pbmcapply::pbmclapply(seq_along(motif_pos), 
                                        function(i){
                                          motif_peaks = rep(0, length(peaks))
                                          motifs = motif_pos[[i]]
                                          overlaps = findOverlaps(motifs, binding_sites_Q17_ranges)
                                          motifs = motifs[overlaps@from]
                                          motif_peaks[unique(findOverlaps(peaks, motifs)@from)] = 1
                                          return(motif_peaks)
                                        }, mc.cores = 2)
motif_fp_matrix <- do.call(cbind, motif_fp_matrix)
motif_fp_matrix = Matrix::Matrix(motif_fp_matrix, sparse = TRUE)
colnames(motif_fp_matrix) = names(pwm)
rownames(motif_fp_matrix) = paste0(regionsBed$chr, '-', 
                                   regionsBed$start, '-', 
                                   regionsBed$end)
weeks56 = seuratobject %>%
  subset(dataset == 'weeks_56')
weeks56 <- RunChromVAR(
  object = weeks56,
  assay = 'ATAC',
  motif.matrix = motif_fp_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks56[['chromvar']], "age_56_weeks/Q17/chromvar_footprints_vertebrates_0.05_120.rds")

weeks56 <- RunChromVAR(
  object = weeks56,
  assay = 'ATAC',
  motif.matrix = motif_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks56[['chromvar']], "age_56_weeks/Q17/chromvar_vertebrates.rds")


### 82 weeks =====
motif_fp_matrix = pbmcapply::pbmclapply(seq_along(motif_pos), 
                                        function(i){
                                          motif_peaks = rep(0, length(peaks))
                                          motifs = motif_pos[[i]]
                                          overlaps = findOverlaps(motifs, binding_sites_U21_ranges)
                                          motifs = motifs[overlaps@from]
                                          motif_peaks[unique(findOverlaps(peaks, motifs)@from)] = 1
                                          return(motif_peaks)
                                        }, 
                                        mc.cores = 2)
motif_fp_matrix <- do.call(cbind, motif_fp_matrix)
motif_fp_matrix = Matrix::Matrix(motif_fp_matrix, sparse = TRUE)
colnames(motif_fp_matrix) = names(pwm)
rownames(motif_fp_matrix) = paste0(regionsBed$chr, '-', 
                                   regionsBed$start, '-', 
                                   regionsBed$end)
weeks82 = seuratobject %>%
  subset(dataset == 'weeks_82')
weeks82 <- RunChromVAR(
  object = weeks82,
  assay = 'ATAC',
  motif.matrix = motif_fp_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks82[['chromvar']], "age_82_weeks/U21/chromvar_footprints_vertebrates_0.05_120.rds")

weeks82 <- RunChromVAR(
  object = weeks82,
  assay = 'ATAC',
  motif.matrix = motif_matrix, 
  genome = genome_file, 
  verbose = TRUE
)

saveRDS(weeks82[['chromvar']], "age_82_weeks/U21/chromvar_vertebrates.rds")


ap1_tf = c("JUN", "JUNB", 
           "JUND", 
           "JUN::JUNB", 
           "FOS", "FOS::JUN", 
           "FOS::JUNB", 
           "FOS::JUND", 
           "FOSB::JUN", 
           "FOSB::JUNB", 
           "FOSL1", "FOSL1::JUN", 
           "FOSL1::JUNB", 
           "BATF::JUN","BATF3", 
           "FOSL1::JUND", 
           "ATF2", 
           "ATF3", "ATF4", "BATF",
           "FOSL2", "FOSL2::JUN", 
           "FOSL2::JUNB", "FOSL2::JUND")
cell_identity_tf = c('HNF1A', 'HNF1B', 
                     'PPARA::RXRA', 
                     'PPARG', 
                     'HNF4G', 'Pparg::Rxra', 
                     'SOX9', 'NR3C1', 'NR3C2', 
                     'NFAT5', 'HNF4A(var.2)', 
                     "TFAP2A", "TFAP2B", 
                     'PAX2',  
                     'RREB', 
                     'HNF4A',
                     "NFKB1", "NFKB2"
)
### 16 .82 comparison ----
ap1_tf = c("JUN", "JUNB", 
           "JUND", 
           "JUN::JUNB", 
           "FOS", "FOS::JUN", 
           "FOS::JUNB", 
           "FOS::JUND", 
           "FOSB::JUN", 
           "FOSB::JUNB", 
           "FOSL1", "FOSL1::JUN", 
           "FOSL1::JUNB", 
           "BATF::JUN","BATF3", 
           "FOSL1::JUND", 
           "ATF2", 
           "ATF3", "ATF4", "BATF",
           "FOSL2", "FOSL2::JUN", 
           "FOSL2::JUNB", "FOSL2::JUND")

weeks82_chromvar_fp = readRDS("age_82_weeks/U21/chromvar_footprints_vertebrates_0.05_120.rds")
weeks16_chromvar_fp = readRDS("age_16_weeks/C03/chromvar_footprints_vertebrates_0.05_120.rds")
weeks82_chromvar = readRDS("age_82_weeks/U21/chromvar_vertebrates.rds")
weeks16_chromvar = readRDS("age_16_weeks/C03/chromvar_vertebrates.rds")

motif_names = names(pwm)
gene_names_all = unlist(sapply(motif_names, function(x) {
  getMatrixByID(JASPAR2020, ID = x)@name
}))
weeks82_PT = rownames(weeks82@meta.data)[weeks82$celltype %in% c('PT', 'PT-MT')]
weeks16_PT = rownames(weeks16@meta.data)[weeks16$celltype %in% c('PT', 'PT-MT')]
## fp chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar_fp@data[tf, weeks82_PT]
  weeks16_value = weeks16_chromvar_fp@data[tf, weeks16_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks16_value))
  U <- wilcox.test(weeks82_value, 
                   weeks16_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks16_value)$p.value
  )
}

wilcox_tf_82_16_fp_chromvar = data.frame('stat' = wilcox_stat, 
                                         'pval' = wilcox_pval, 
                                         'tf' = gene_names_all)
wilcox_tf_82_16_fp_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_16_fp_chromvar$pval, 'fdr')
wilcox_tf_82_16_fp_chromvar[wilcox_tf_82_16_fp_chromvar$tf %in% ap1_tf, ]

## chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar@data[tf, weeks82_PT]
  weeks16_value = weeks16_chromvar@data[tf, weeks16_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks16_value))
  U <- wilcox.test(weeks82_value, 
                   weeks16_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks16_value)$p.value
  )
}

wilcox_tf_82_16_chromvar = data.frame('stat' = wilcox_stat, 
                                      'pval' = wilcox_pval, 
                                      'tf' = gene_names_all)
wilcox_tf_82_16_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_16_chromvar$pval, 'fdr')
wilcox_tf_82_16_chromvar[wilcox_tf_82_16_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_82_16_chromvar[wilcox_tf_82_16_chromvar$tf %in% ap1_tf, ]
### diff plot -----
wilcox_tf_82_16_fp_chromvar = wilcox_tf_82_16_fp_chromvar[order(wilcox_tf_82_16_fp_chromvar$stat, 
                                                                decreasing = F), ]
wilcox_tf_82_16_fp_chromvar$rank = 1:nrow(wilcox_tf_82_16_fp_chromvar)
write.csv(wilcox_tf_82_16_fp_chromvar, 
          '16v.s.82/differential_chromvar_fp_120_0.05.csv')

wilcox_tf_82_16_chromvar = wilcox_tf_82_16_chromvar[order(wilcox_tf_82_16_chromvar$stat, 
                                                          decreasing = F), ]
wilcox_tf_82_16_chromvar$rank = 1:nrow(wilcox_tf_82_16_chromvar)
write.csv(wilcox_tf_82_16_chromvar, 
          '16v.s.82/differential_chromvar.csv')

gene_interest = c('CTCF', "NRF1", 'NFKB1', 
                  # 'STAT1', # 'Stat5a', # 'KLF9', 
                  # 'TEAD1', 
                  # 'Stat4', 
                  # 'STAT3', # 'Stat5b', 
                  # 'NR3C1', 
                  # 'KLF15', 
                  'BATF3', 'FOSL1', 
                  'BATF', 'FOS', 'FOSL2', 
                  'JUN(var.2)', 'JUND', 'JUNB', 
                  'HNF4A', # 'Stat2', 
                  'IRF1', # 'Stat5a::Stat5b', 
                  # 'Klf12', # 'FOSB::JUN', 
                  # 'JUN', 
                  'ATF4'# , # 'ATF3', 
                  # 'FOS::JUNB'# , 'STAT1::STAT2'
)
wilcox_tf_82_16_fp_chromvar = wilcox_tf_82_16_fp_chromvar %>%
  mutate('test_result' = ifelse(pval.adjusted <= 0.05, 
                                ifelse(stat < 0, 'down regulation', 
                                       'up regulation'), 
                                'insignificant'))
p1 = ggplot(wilcox_tf_82_16_fp_chromvar, aes(x = rank, 
                                             y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', 
             linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_16_fp_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_16_fp_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 21, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar combined footprint") + 
  scale_x_continuous(labels = NULL, breaks = NULL,
                     limits = c(0, 747)) + 
  labs(x = '') +
  theme(legend.title=element_blank())
p1
ggsave('16v.s.82/chromvar_fp_rank_120_0.05.pdf', p1, 
       height = 3, width = 4)

ggplot(wilcox_tf_82_16_chromvar, aes(x = rank, 
                                     y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_16_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_16_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(1, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 30, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar") + 
  scale_x_continuous(labels = NULL, breaks = NULL, 
                     limits = c(0, 747)) + 
  xlab('')

ggsave('16v.s.82/chromvar_rank.pdf', 
       height = 3, width = 4)

### 30 .56 comparison ----
ap1_tf = c("JUN", "JUNB", 
           "JUND", 
           "JUN::JUNB", 
           "FOS", "FOS::JUN", 
           "FOS::JUNB", 
           "FOS::JUND", 
           "FOSB::JUN", 
           "FOSB::JUNB", 
           "FOSL1", "FOSL1::JUN", 
           "FOSL1::JUNB", 
           "BATF::JUN","BATF3", 
           "FOSL1::JUND", 
           "ATF2", 
           "ATF3", "ATF4", "BATF",
           "FOSL2", "FOSL2::JUN", 
           "FOSL2::JUNB", "FOSL2::JUND")

weeks56_chromvar_fp = readRDS("age_56_weeks/Q17/chromvar_footprints_vertebrates_0.05_120.rds")
weeks30_chromvar_fp = readRDS("age_30_weeks/L12/chromvar_footprints_vertebrates_0.05_120.rds")
weeks56_chromvar = readRDS("age_56_weeks/Q17/chromvar_vertebrates.rds")
weeks30_chromvar = readRDS("age_30_weeks/L12/chromvar_vertebrates.rds")

motif_names = names(pwm)
gene_names_all = unlist(sapply(motif_names, function(x) {
  getMatrixByID(JASPAR2020, ID = x)@name
}))
weeks56_PT = rownames(weeks56@meta.data)[weeks56$celltype %in% c('PT', 'PT-MT')]
weeks30_PT = rownames(weeks30@meta.data)[weeks30$celltype %in% c('PT', 'PT-MT')]

## fp chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks56_value = weeks56_chromvar_fp@data[tf, weeks56_PT]
  weeks30_value = weeks30_chromvar_fp@data[tf, weeks30_PT]
  n1 <- sum(!is.na(weeks56_value))
  n2 <- sum(!is.na(weeks30_value))
  U <- wilcox.test(weeks56_value, 
                   weeks30_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks56_value, 
                              weeks30_value)$p.value
  )
}

wilcox_tf_56_30_fp_chromvar = data.frame('stat' = wilcox_stat, 
                                         'pval' = wilcox_pval, 
                                         'tf' = gene_names_all)
wilcox_tf_56_30_fp_chromvar$pval.adjusted = p.adjust(wilcox_tf_56_30_fp_chromvar$pval, 'fdr')
wilcox_tf_56_30_fp_chromvar[wilcox_tf_56_30_fp_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_56_30_fp_chromvar[wilcox_tf_56_30_fp_chromvar$tf %in% ap1_tf, ]

## chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks56_value = weeks56_chromvar@data[tf, weeks56_PT]
  weeks30_value = weeks30_chromvar@data[tf, weeks30_PT]
  n1 <- sum(!is.na(weeks56_value))
  n2 <- sum(!is.na(weeks30_value))
  U <- wilcox.test(weeks56_value, 
                   weeks30_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks56_value, 
                              weeks30_value)$p.value
  )
}

wilcox_tf_56_30_chromvar = data.frame('stat' = wilcox_stat, 
                                      'pval' = wilcox_pval, 
                                      'tf' = gene_names_all)
wilcox_tf_56_30_chromvar$pval.adjusted = p.adjust(wilcox_tf_56_30_chromvar$pval, 'fdr')
wilcox_tf_56_30_chromvar[wilcox_tf_56_30_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_56_30_chromvar[wilcox_tf_56_30_chromvar$tf %in% ap1_tf, ]

### diff plot -----
wilcox_tf_56_30_fp_chromvar = wilcox_tf_56_30_fp_chromvar[order(wilcox_tf_56_30_fp_chromvar$stat, 
                                                                decreasing = F), ]
wilcox_tf_56_30_fp_chromvar$rank = 1:nrow(wilcox_tf_56_30_fp_chromvar)
write.csv(wilcox_tf_56_30_fp_chromvar, 
          '30v.s.56/differential_chromvar_fp_120.csv')

wilcox_tf_56_30_chromvar = wilcox_tf_56_30_chromvar[order(wilcox_tf_56_30_chromvar$stat, 
                                                          decreasing = F), ]
wilcox_tf_56_30_chromvar$rank = 1:nrow(wilcox_tf_56_30_chromvar)
write.csv(wilcox_tf_56_30_chromvar, 
          '30v.s.56/differential_chromvar.csv')

gene_interest = c('CTCF', "NRF1", 'NFKB1', 
                  'STAT1', # 'Stat5a', # 'KLF9', 
                  # 'TEAD1', 
                  # 'Stat4', 
                  # 'STAT3', # 'Stat5b', 
                  'NR3C1', 
                  # 'KLF15', 
                  'BATF3', 'FOS::JUN', 
                  'FOSB::JUNB', 'FOSL1::JUNB', 
                  'BATF', 'FOS', 'FOSL2', 
                  'JUN(var.2)', # 'JUND', 'JUNB', 
                  'HNF4A', # 'Stat2', 
                  # 'IRF1', # 'Stat5a::Stat5b', 
                  # 'Klf12', # 'FOSB::JUN', 
                  # 'JUN', 
                  'ATF4'# , # 'ATF3', 
                  # 'FOS::JUNB'# , 'STAT1::STAT2'
)
wilcox_tf_56_30_fp_chromvar = wilcox_tf_56_30_fp_chromvar %>%
  mutate('test_result' = ifelse(pval.adjusted <= 0.05, 
                                ifelse(stat < 0, 'down regulation', 
                                       'up regulation'), 
                                'insignificant'))
p1 = ggplot(wilcox_tf_56_30_fp_chromvar, aes(x = rank, 
                                             y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', 
             linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_56_30_fp_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_56_30_fp_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 20, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar combined footprint") + 
  scale_x_continuous(labels = NULL, breaks = NULL,
                     limits = c(0, 747)) + 
  labs(x = '') +
  theme(legend.title=element_blank())
p1
ggsave('30v.s.56/chromvar_fp_rank_120_0.05.pdf', p1, 
       height = 3, width = 4)

ggplot(wilcox_tf_56_30_chromvar, aes(x = rank, 
                                     y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_56_30_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_56_30_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(1, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 40, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar") + 
  scale_x_continuous(labels = NULL, breaks = NULL, 
                     limits = c(0, 747)) + 
  xlab('')

ggsave('30v.s.56/chromvar_rank.pdf', 
       height = 3, width = 4)

### 56 .82 comparison ----
ap1_tf = c("JUN", "JUNB", 
           "JUND", 
           "JUN::JUNB", 
           "FOS", "FOS::JUN", 
           "FOS::JUNB", 
           "FOS::JUND", 
           "FOSB::JUN", 
           "FOSB::JUNB", 
           "FOSL1", "FOSL1::JUN", 
           "FOSL1::JUNB", 
           "BATF::JUN","BATF3", 
           "FOSL1::JUND", 
           "ATF2", 
           "ATF3", "ATF4", "BATF",
           "FOSL2", "FOSL2::JUN", 
           "FOSL2::JUNB", "FOSL2::JUND")

weeks82_chromvar_fp = readRDS("age_82_weeks/U21/chromvar_footprints_vertebrates_0.05_120.rds")
weeks56_chromvar_fp = readRDS("age_56_weeks/Q17/chromvar_footprints_vertebrates_0.05_120.rds")
weeks82_chromvar = readRDS("age_82_weeks/U21/chromvar_vertebrates.rds")
weeks56_chromvar = readRDS("age_56_weeks/Q17/chromvar_vertebrates.rds")

motif_names = names(pwm)
gene_names_all = unlist(sapply(motif_names, function(x) {
  getMatrixByID(JASPAR2020, ID = x)@name
}))
weeks82_PT = rownames(weeks82@meta.data)[weeks82$celltype %in% c('PT', 'PT-MT')]
weeks56_PT = rownames(weeks56@meta.data)[weeks56$celltype %in% c('PT', 'PT-MT')]

## fp chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar_fp@data[tf, weeks82_PT]
  weeks56_value = weeks56_chromvar_fp@data[tf, weeks56_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks56_value))
  U <- wilcox.test(weeks82_value, 
                   weeks56_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks56_value)$p.value
  )
}

wilcox_tf_82_56_fp_chromvar = data.frame('stat' = wilcox_stat, 
                                         'pval' = wilcox_pval, 
                                         'tf' = gene_names_all)
wilcox_tf_82_56_fp_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_56_fp_chromvar$pval, 'fdr')
wilcox_tf_82_56_fp_chromvar[wilcox_tf_82_56_fp_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_82_56_fp_chromvar[wilcox_tf_82_56_fp_chromvar$tf %in% ap1_tf, ]

## chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar@data[tf, weeks82_PT]
  weeks56_value = weeks56_chromvar@data[tf, weeks56_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks56_value))
  U <- wilcox.test(weeks82_value, 
                   weeks56_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks56_value)$p.value
  )
}

wilcox_tf_82_56_chromvar = data.frame('stat' = wilcox_stat, 
                                      'pval' = wilcox_pval, 
                                      'tf' = gene_names_all)
wilcox_tf_82_56_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_56_chromvar$pval, 'fdr')
wilcox_tf_82_56_chromvar[wilcox_tf_82_56_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_82_56_chromvar[wilcox_tf_82_56_chromvar$tf %in% ap1_tf, ]
### diff plot -----
wilcox_tf_82_56_fp_chromvar = wilcox_tf_82_56_fp_chromvar[order(wilcox_tf_82_56_fp_chromvar$stat, 
                                                                decreasing = F), ]
wilcox_tf_82_56_fp_chromvar$rank = 1:nrow(wilcox_tf_82_56_fp_chromvar)
write.csv(wilcox_tf_82_56_fp_chromvar, 
          '56v.s.82/differential_chromvar_fp_120.csv')

wilcox_tf_82_56_chromvar = wilcox_tf_82_56_chromvar[order(wilcox_tf_82_56_chromvar$stat, 
                                                          decreasing = F), ]
wilcox_tf_82_56_chromvar$rank = 1:nrow(wilcox_tf_82_56_chromvar)
write.csv(wilcox_tf_82_56_chromvar, 
          '56v.s.82/differential_chromvar.csv')

gene_interest = c('CTCF', "NRF1", 'NFKB1', 
                  'STAT1', # 'Stat5a', # 'KLF9', 
                  # 'TEAD1', 
                  # 'Stat4', 
                  # 'STAT3', # 'Stat5b', 
                  'NR3C1', 
                  # 'KLF15', 
                  'BATF', 'FOS', 'FOSL2', 
                  'JUN(var.2)', 'JUND', 'JUNB', 
                  'HNF4A', 'Stat2', 
                  # 'IRF1', # 'Stat5a::Stat5b', 
                  # 'Klf12', # 'FOSB::JUN', 
                  # 'JUN', 
                  'ATF4'# , # 'ATF3', 
                  # 'FOS::JUNB'# , 'STAT1::STAT2'
)
wilcox_tf_82_56_fp_chromvar = wilcox_tf_82_56_fp_chromvar %>%
  mutate('test_result' = ifelse(pval.adjusted <= 0.05, 
                                ifelse(stat < 0, 'down regulation', 
                                       'up regulation'), 
                                'insignificant'))
p1 = ggplot(wilcox_tf_82_56_fp_chromvar, aes(x = rank, 
                                             y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', 
             linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_56_fp_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_56_fp_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 20, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar combined footprint") + 
  scale_x_continuous(labels = NULL, breaks = NULL,
                     limits = c(0, 747)) + 
  labs(x = '') +
  theme(legend.title=element_blank())
p1
ggsave('56v.s.82/chromvar_fp_rank_120_0.05.pdf', p1, 
       height = 3, width = 4)

ggplot(wilcox_tf_82_56_chromvar, aes(x = rank, 
                                     y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_56_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_56_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(1, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 25, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar") + 
  scale_x_continuous(labels = NULL, breaks = NULL, 
                     limits = c(0, 747)) + 
  xlab('')

ggsave('56v.s.82/chromvar_rank.pdf', 
       height = 3, width = 4)

### 30 .82 comparison ----
ap1_tf = c("JUN", "JUNB", 
           "JUND", 
           "JUN::JUNB", 
           "FOS", "FOS::JUN", 
           "FOS::JUNB", 
           "FOS::JUND", 
           "FOSB::JUN", 
           "FOSB::JUNB", 
           "FOSL1", "FOSL1::JUN", 
           "FOSL1::JUNB", 
           "BATF::JUN","BATF3", 
           "FOSL1::JUND", 
           "ATF2", 
           "ATF3", "ATF4", "BATF",
           "FOSL2", "FOSL2::JUN", 
           "FOSL2::JUNB", "FOSL2::JUND")

weeks82_chromvar_fp = readRDS("age_82_weeks/U21/chromvar_footprints_vertebrates_0.05_120.rds")
weeks30_chromvar_fp = readRDS("age_30_weeks/L12/chromvar_footprints_vertebrates_0.05_120.rds")
weeks82_chromvar = readRDS("age_82_weeks/U21/chromvar.rds")
weeks30_chromvar = readRDS("age_30_weeks/L12/chromvar.rds")

motif_names = names(pwm)
gene_names_all = unlist(sapply(motif_names, function(x) {
  getMatrixByID(JASPAR2020, ID = x)@name
}))
weeks82_PT = rownames(weeks82@meta.data)[weeks82$celltype %in% c('PT', 'PT-MT')]
weeks30_PT = rownames(weeks30@meta.data)[weeks30$celltype %in% c('PT', 'PT-MT')]

## fp chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar_fp@data[tf, weeks82_PT]
  weeks30_value = weeks30_chromvar_fp@data[tf, weeks30_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks30_value))
  U <- wilcox.test(weeks82_value, 
                   weeks30_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks30_value)$p.value
             )
}

wilcox_tf_82_30_fp_chromvar = data.frame('stat' = wilcox_stat, 
                                         'pval' = wilcox_pval, 
                                         'tf' = gene_names_all)
wilcox_tf_82_30_fp_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_30_fp_chromvar$pval, 'fdr')
wilcox_tf_82_30_fp_chromvar[wilcox_tf_82_30_fp_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_82_30_fp_chromvar[wilcox_tf_82_30_fp_chromvar$tf %in% ap1_tf, ]

## chromvar ----
wilcox_stat = NULL
wilcox_pval = NULL
for (tf in motif_names){
  weeks82_value = weeks82_chromvar@data[tf, weeks82_PT]
  weeks30_value = weeks30_chromvar@data[tf, weeks30_PT]
  n1 <- sum(!is.na(weeks82_value))
  n2 <- sum(!is.na(weeks30_value))
  U <- wilcox.test(weeks82_value, 
                   weeks30_value)$statistic
  mu <- n1 * n2 / 2
  sigma <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  z <- (U - mu) / sigma
  
  wilcox_stat = c(wilcox_stat, 
                  z)
  wilcox_pval = c(wilcox_pval,
                  wilcox.test(weeks82_value, 
                              weeks30_value)$p.value
  )
}

wilcox_tf_82_30_chromvar = data.frame('stat' = wilcox_stat, 
                                      'pval' = wilcox_pval, 
                                      'tf' = gene_names_all)
wilcox_tf_82_30_chromvar$pval.adjusted = p.adjust(wilcox_tf_82_30_chromvar$pval, 'fdr')
wilcox_tf_82_30_chromvar[wilcox_tf_82_30_chromvar$tf %in% cell_identity_tf, ]
wilcox_tf_82_30_chromvar[wilcox_tf_82_30_chromvar$tf %in% ap1_tf, ]

### diff plot -----
wilcox_tf_82_30_fp_chromvar = wilcox_tf_82_30_fp_chromvar[order(wilcox_tf_82_30_fp_chromvar$stat, 
                                                                decreasing = F), ]
wilcox_tf_82_30_fp_chromvar$rank = 1:nrow(wilcox_tf_82_30_fp_chromvar)
write.csv(wilcox_tf_82_30_fp_chromvar, 
          '30v.s.82/differential_chromvar_fp_120.csv')

wilcox_tf_82_30_chromvar = wilcox_tf_82_30_chromvar[order(wilcox_tf_82_30_chromvar$stat, 
                                                                decreasing = F), ]
wilcox_tf_82_30_chromvar$rank = 1:nrow(wilcox_tf_82_30_chromvar)
write.csv(wilcox_tf_82_30_chromvar, 
          '30v.s.82/differential_chromvar.csv')

gene_interest = c('CTCF', "NRF1", 'NFKB1', 
                  'STAT1', # 'Stat5a', # 'KLF9', 
                  # 'TEAD1', 
                  # 'Stat4', 
                  # 'STAT3', # 'Stat5b', 
                  'NR3C1', 
                  # 'KLF15', 
                  'BATF', 'FOS', 'FOSL2', 
                  'JUN(var.2)', 'JUND', 'JUNB', 
                  'HNF4A', 'Stat2', 
                  'IRF1', # 'Stat5a::Stat5b', 
                  # 'Klf12', # 'FOSB::JUN', 
                  # 'JUN', 
                  'ATF4'# , # 'ATF3', 
                  # 'FOS::JUNB'# , 'STAT1::STAT2'
                  )
wilcox_tf_82_30_fp_chromvar = wilcox_tf_82_30_fp_chromvar %>%
  mutate('test_result' = ifelse(pval.adjusted <= 0.05, 
                                ifelse(stat < 0, 'down regulation', 
                                       'up regulation'), 
                                'insignificant'))
p1 = ggplot(wilcox_tf_82_30_fp_chromvar, aes(x = rank, 
                                        y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', 
             linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_30_fp_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_30_fp_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 20, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar combined footprint") + 
  scale_x_continuous(labels = NULL, breaks = NULL,
                     limits = c(0, 747)) + 
  labs(x = '') +
  theme(legend.title=element_blank())
p1
ggsave('30v.s.82/chromvar_fp_rank_120_0.05.pdf', p1, 
       height = 3, width = 4)

ggplot(wilcox_tf_82_30_chromvar, aes(x = rank, 
                                        y = stat)) + 
  geom_hline(yintercept = 0, color = 'gray', linetype = 'dashed') +
  geom_point(size = 0.6, color = '#9ED2BE', alpha = 0.3) + 
  geom_point(data = subset(wilcox_tf_82_30_chromvar, 
                           tf %in% gene_interest), 
             size = 0.8, color = 'red', alpha = 0.9) + 
  ylab("normalized Wilcoxon stat") + 
  ggrepel::geom_text_repel(data = subset(wilcox_tf_82_30_chromvar, 
                                         tf %in% gene_interest), 
                           color = 'black',
                           aes(label = tf), 
                           box.padding = unit(1, 'lines'),
                           point.padding = unit(0.3, 'lines'),
                           max.overlaps = 25, size = 3.5) +
  theme_Publication() + 
  ggtitle("chromVar") + 
  scale_x_continuous(labels = NULL, breaks = NULL, 
                     limits = c(0, 747)) + 
  xlab('')

ggsave('30v.s.82/chromvar_rank.pdf', 
       height = 3, width = 4)



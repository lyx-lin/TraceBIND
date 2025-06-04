library(tidyverse)
library(GenomicRanges)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all"))
source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')
### 16 weeks -----
files_16 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_16_weeks/C03/effect_size", 
                      full.names = FALSE)
files_30 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_30_weeks/L12/effect_size", 
                      full.names = FALSE)
files_56 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_56_weeks/Q17/effect_size", 
                      full.names = FALSE)
files_82 = list.files("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_82_weeks/U21/effect_size", 
                      full.names = FALSE)
files = Reduce(intersect, list(files_16,
                               files_30,
                               files_56, 
                               files_82
))
files = gtools::mixedsort(files)

logp_u21_multiome = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/logp_threshold.txt"))
logp_q17_multiome = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_56_weeks/Q17/kidney_multiome/logp_threshold.txt"))
logp_l12_multiome = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_30_weeks/L12/kidney_multiome/logp_threshold.txt"))
logp_c03_multiome = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_16_weeks/C03/kidney_multiome/logp_threshold.txt"))

setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_16_weeks/C03/effect_size")
binding_sites_C03 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_c03_multiome$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(rownames(logp_c03_multiome))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)
# binding_sites_C03 = data.table::rbindlist(binding_sites_C03)
saveRDS(binding_sites_C03, 
        '/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/age_16_weeks/C03/thresholded_binding_sites.rds')

### 30 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_30_weeks/L12/effect_size")
binding_sites_L12 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_l12_multiome$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                     as.numeric(rownames(logp_l12_multiome))))]
                                              )
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)
# binding_sites_L12 = data.table::rbindlist(binding_sites_L12)
saveRDS(binding_sites_L12, 
        '/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/age_30_weeks/L12/thresholded_binding_sites_total.rds')

### 56 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_56_weeks/Q17/effect_size")
binding_sites_Q17 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$total_coverage = binding_sites$mean_coverage*binding_sites$width
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_q17_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_q17_multiome_mean$labels_mean_counts))))]
                                                )
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)
# binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17)
saveRDS(binding_sites_Q17, 
        '/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/age_56_weeks/Q17/thresholded_binding_sites_downsample.rds')

### 82 weeks ----
setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/age_82_weeks/U21/effect_size")
binding_sites_U21 = pbmcapply::pbmclapply(1:length(files), 
                                          function(i) {
                                            binding_sites = read.table(files[i], header = TRUE) 
                                            j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                            if (dim(binding_sites)[1] > 0){
                                              binding_sites$total_coverage = binding_sites$mean_coverage*binding_sites$width
                                              binding_sites$threshold = apply(binding_sites, 1, function(x)
                                                logp_u21_multiome_mean$smoothed_threshold[which.min(abs(x['mean_coverage'] - 
                                                                                                          as.numeric(as.character(logp_u21_multiome_mean$labels_mean_counts))))])
                                              binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                              binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                            }
                                            binding_sites
                                          }, 
                                          mc.cores = 2
)

# binding_sites_U21 = data.table::rbindlist(binding_sites_U21)
saveRDS(binding_sites_U21, 
        '/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/age_82_weeks/U21/thresholded_binding_sites.rds')

### visualizations ----
setwd("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")
binding_sites_C03 = readRDS('age_16_weeks/C03/thresholded_binding_sites.rds')
binding_sites_L12 = readRDS('age_30_weeks/L12/thresholded_binding_sites.rds')
binding_sites_Q17 = readRDS('age_56_weeks/Q17/thresholded_binding_sites.rds')
binding_sites_U21 = readRDS('age_82_weeks/U21/thresholded_binding_sites.rds')

binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill=TRUE)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill=TRUE)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill=TRUE)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill=TRUE)

count_footprinting_C03 = NULL
count_footprinting_L12 = NULL
count_footprinting_Q17 = NULL
count_footprinting_U21 = NULL

widths = seq(15, 200, 10)
for (j in widths){
  count_footprinting_C03 = c(count_footprinting_C03, 
                             dim(binding_sites_C03[binding_sites_C03$width_effect_size < j + 8 & 
                                                     binding_sites_C03$width_effect_size >= j - 6, ])[1])
  count_footprinting_L12 = c(count_footprinting_L12, 
                             dim(binding_sites_L12[binding_sites_L12$width_effect_size < j + 8 & 
                                                     binding_sites_L12$width_effect_size >= j - 6])[1])
  count_footprinting_Q17 = c(count_footprinting_Q17, 
                             dim(binding_sites_Q17[binding_sites_Q17$width_effect_size < j + 8 & 
                                                     binding_sites_Q17$width_effect_size >= j - 6])[1])
  count_footprinting_U21 = c(count_footprinting_U21, 
                             dim(binding_sites_U21[binding_sites_U21$width_effect_size < j + 8 & 
                                                     binding_sites_U21$width_effect_size >= j - 6])[1])
  
}
peaks_width = regionsBed$end - regionsBed$start + 1

# widths = c(15, widths)
results = data.frame('width' = widths, 
                              'weeks_16' = count_footprinting_C03*1000/sum(peaks_width), 
                              'weeks_30' = count_footprinting_L12*1000/sum(peaks_width), 
                              'weeks_56' = count_footprinting_Q17*1000/sum(peaks_width), 
                              'weeks_82' = count_footprinting_U21*1000/sum(peaks_width)
                     )
results = pivot_longer(results, 
                       cols = c('weeks_16', 
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82'
                                ),  
                       names_to = 'age',         
                       values_to = 'count')
results$age = factor(results$age, 
                     levels = c('weeks_16',
                                'weeks_30', 
                                'weeks_56', 
                                'weeks_82')
                     )

p1 = ggplot(results) + 
  geom_point(aes(x = width, y = count, color = age), 
             size = 1.5) + 
  geom_line(aes(x = width, y = count, color = age), 
            linewidth = 1.5) + 
  ylab("footprintings number per kb") + 
  scale_x_continuous(limits = c(0, 200), 
                     breaks = c(0, seq(50, 200, 50))) + 
  labs(colour = "Age") + 
  scale_color_manual(values = c('#D3E671',
                                '#A3D8FF', 
                                '#FAA301',
                                '#E58C94'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks"))
p1
ggsave("age_num_comparison_subset.png", p1, 
       height = 4, width = 6)

### merging seuratobject ----
library(Signac)
library(Seurat)
metadata = read.table("~/nzhanglab/data/ParkerWilson_DEFND/metadata.txt")
# metadata = metadata[metadata$library_type == 'Multiome',]
cells_conversion = read.csv("~/nzhanglab/data/ParkerWilson_DEFND/barcode_conversion_10Xmultiome.csv")
peaks = makeGRangesFromDataFrame(regionsBed)
### 16 weeks ----
metadata_16 = metadata[metadata$library_id == "RAGE24-C03-KYC-LN-01-Multiome-RNA", ]
metadata_16 = metadata_16['celltype']
metadata_16$GEX_bc = rownames(metadata_16)
metadata_16$GEX_bc = gsub(substr(metadata_16$GEX_bc[1], nchar(metadata_16$GEX_bc[1]) + 1 - 2, nchar(metadata_16$GEX_bc[1])), '', metadata_16$GEX_bc)

barcodeGroups_16 = left_join(metadata_16, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_16$GEX_bc = NULL
barcodeGroups_16 = barcodeGroups_16[c(2, 1)]
barcodeGroups_16$ATAC_bc = paste0(barcodeGroups_16$ATAC_bc, '-1')
frags.16 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-C03-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_16$ATAC_bc
)

### 30 weeks -----
metadata_30 = metadata[metadata$library_id == "RAGE24-L12-KYC-LN-01-Multiome-RNA", ]
metadata_30 = metadata_30['celltype']
metadata_30$GEX_bc = rownames(metadata_30)
metadata_30$GEX_bc = gsub(substr(metadata_30$GEX_bc[1], nchar(metadata_30$GEX_bc[1]) + 1 - 3, nchar(metadata_30$GEX_bc[1])), '', metadata_30$GEX_bc)

barcodeGroups_30 = left_join(metadata_30, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_30$GEX_bc = NULL
barcodeGroups_30 = barcodeGroups_30[c(2, 1)]
barcodeGroups_30$ATAC_bc = paste0(barcodeGroups_30$ATAC_bc, '-1')
frags.30 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-L12-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_30$ATAC_bc
)

### 56 weeks -----
metadata_56 = metadata[metadata$library_id == "RAGE24-Q17-KYC-LN-01-Multiome-RNA", ]
metadata_56 = metadata_56['celltype']
metadata_56$GEX_bc = rownames(metadata_56)
metadata_56$GEX_bc = gsub(substr(metadata_56$GEX_bc[1], nchar(metadata_56$GEX_bc[1]) + 1 - 3, nchar(metadata_56$GEX_bc[1])), '', metadata_56$GEX_bc)

barcodeGroups_56 = left_join(metadata_56, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_56$GEX_bc = NULL
barcodeGroups_56 = barcodeGroups_56[c(2, 1)]
barcodeGroups_56$ATAC_bc = paste0(barcodeGroups_56$ATAC_bc, '-1')
frags.56 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-Q17-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_56$ATAC_bc
)

### 82 weeks ----
metadata_82 = metadata[metadata$library_id == "RAGE24-U21-KYC-LN-01-Multiome-RNA", ]
metadata_82 = metadata_82['celltype']
metadata_82$GEX_bc = rownames(metadata_82)
metadata_82$GEX_bc = gsub(substr(metadata_82$GEX_bc[1], nchar(metadata_82$GEX_bc[1]) + 1 - 3, nchar(metadata_82$GEX_bc[1])), '', metadata_82$GEX_bc)

barcodeGroups_82 = left_join(metadata_82, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_82$GEX_bc = NULL
barcodeGroups_82 = barcodeGroups_82[c(2, 1)]
barcodeGroups_82$ATAC_bc = paste0(barcodeGroups_82$ATAC_bc, '-1')
frags.82 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_82$ATAC_bc
)
### combining -----
weeks16.counts <- FeatureMatrix(
  fragments = frags.16,
  features = peaks,
  cells = barcodeGroups_16$ATAC_bc
)

weeks30.counts <- FeatureMatrix(
  fragments = frags.30,
  features = peaks,
  cells = barcodeGroups_30$ATAC_bc
)

weeks56.counts <- FeatureMatrix(
  fragments = frags.56,
  features = peaks,
  cells = barcodeGroups_56$ATAC_bc
)

weeks82.counts <- FeatureMatrix(
  fragments = frags.82,
  features = peaks,
  cells = barcodeGroups_82$ATAC_bc
)

weeks16_assay <- CreateChromatinAssay(weeks16.counts, 
                                      fragments = frags.16)
weeks16 <- CreateSeuratObject(weeks16_assay, 
                              assay = "ATAC", 
                              meta.data=barcodeGroups_16)

weeks30_assay <- CreateChromatinAssay(weeks30.counts, 
                                      fragments = frags.30)
weeks30 <- CreateSeuratObject(weeks30_assay, 
                              assay = "ATAC", 
                              meta.data=barcodeGroups_30)

weeks56_assay <- CreateChromatinAssay(weeks56.counts, 
                                      fragments = frags.56)
weeks56 <- CreateSeuratObject(weeks56_assay, 
                              assay = "ATAC", 
                              meta.data=barcodeGroups_56)

weeks82_assay <- CreateChromatinAssay(weeks82.counts, 
                                      fragments = frags.82)
weeks82 <- CreateSeuratObject(weeks82_assay, 
                              assay = "ATAC", 
                              meta.data=barcodeGroups_82)

weeks16$dataset <- 'weeks_16'
weeks30$dataset <- 'weeks_30'
weeks56$dataset <- 'weeks_56'
weeks82$dataset <- 'weeks_82'

combined <- merge(
  x = weeks16,
  y = list(weeks30, weeks56, weeks82),
  add.cell.ids = c("weeks_16", "weeks_30", "weeks_56", "weeks_82")
)
  
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
DimPlot(combined, group.by = 'celltype', pt.size = 0.1)
saveRDS(combined, 
        'seuratobject_whole_peaks.rds')

seuratobject = readRDS("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/seuratobject_whole_peaks.rds")
seuratobject <- RunTFIDF(seuratobject)
seuratobject <- FindTopFeatures(seuratobject, 
                                min.cutoff = 20)
seuratobject <- RunSVD(seuratobject)

seuratobject <- harmony::RunHarmony(
  object = seuratobject,
  group.by.vars = 'dataset',
  assay.use = 'ATAC',
  reduction.use = 'lsi',
  reduction.save = 'harmony',
  project.dim = FALSE
)

seuratobject <- RunUMAP(seuratobject, 
                        dims = 2:30, 
                        reduction = 'harmony', 
                        reduction.name = 'umap_harmony'
)

beach = c("pink", '#99CC66', "#00A08A", "#0B775E", 
          "#87D2DB","#5BB1CB",'#88CEEF', '#1E78B4', "#4F66AF", "#352A86", 
          '#a8ddb5',"#FF6800","#FCBF6E", "#fbdf72", 
          "#8c6bb1", "#7F3F98", "#FCB31A", '#E6C2DC')
p2 = ggplot() + 
  geom_point(aes(x = Embeddings(seuratobject[['umap_harmony']])[, 1], 
                 y = Embeddings(seuratobject[['umap_harmony']])[, 2], 
                 color = seuratobject$celltype), 
             size = 0.8) + 
  scale_colour_manual(values = beach, 
                      name = 'cell type') + 
  ylab("umap_2") + 
  xlab("umap_1") + 
  theme_Publication(base_size = 14)
p2
ggsave("umap_integration_clustering_whole_peaks.pdf", 
       p2, 
       height = 5, width = 7)
summerNight = c("#FCBF6E", "#bed678","#CAB2D6","#A1CDE1")

p2 = ggplot() + 
  geom_point(aes(x = Embeddings(seuratobject[['umap_harmony']])[, 1], 
                 y = Embeddings(seuratobject[['umap_harmony']])[, 2], 
                 color = seuratobject$dataset), 
             size = 1) +  
  scale_colour_manual(values = summerNight, 
                      name = 'ages') + 
  ylab("umap_2") + 
  xlab("umap_1") + 
  theme_Publication()
ggsave("umap_integration_clustering_whole_peaks_samples.pdf", 
       p2, 
       height = 5, width = 7)

### loci -----
setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all"))
mouse_clock_by_human <- readRDS('../mouse_clock_lifted_human_to_mouse_mm10.rds')
rat_clock_by_human <- easyLift::easyLiftOver(mouse_clock_by_human, map = '../mm10ToRn7.over.chain')
rat_clock_by_human = data.frame(rat_clock_by_human)
reference = read.csv("~/nzhanglab/data/ParkerWilson_DEFND/refernece/GRCr8_GCF_036323735.1_sequence_info.csv")
reference = reference[c('refseqAccession', 'sequenceName')]
rat_clock_by_human = left_join(rat_clock_by_human, reference, 
                               by = join_by(seqnames == sequenceName))
rat_clock_by_human = rat_clock_by_human %>%
  dplyr::select(refseqAccession, start, end)
colnames(rat_clock_by_human)[1] = 'chr'
rat_clock_by_human = makeGRangesFromDataFrame(rat_clock_by_human)

regionsBed = read.table("all_peaks_region.bed")
peaks = makeGRangesFromDataFrame(regionsBed)
aging_loci = unique(findOverlapPairs(rat_clock_by_human, 
                                     peaks)@first)

binding_sites_C03 = data.table::rbindlist(binding_sites_C03, fill=TRUE)
binding_sites_L12 = data.table::rbindlist(binding_sites_L12, fill=TRUE)
binding_sites_Q17 = data.table::rbindlist(binding_sites_Q17, fill=TRUE)
binding_sites_U21 = data.table::rbindlist(binding_sites_U21, fill=TRUE)

binding_sites_C03$start = binding_sites_C03$position_effect_size - binding_sites_C03$width_effect_size/2
binding_sites_C03$end = binding_sites_C03$position_effect_size + binding_sites_C03$width_effect_size/2
binding_sites_L12$start = binding_sites_L12$position_effect_size - binding_sites_L12$width_effect_size/2
binding_sites_L12$end = binding_sites_L12$position_effect_size + binding_sites_L12$width_effect_size/2
binding_sites_Q17$start = binding_sites_Q17$position_effect_size - binding_sites_Q17$width_effect_size/2
binding_sites_Q17$end = binding_sites_Q17$position_effect_size + binding_sites_Q17$width_effect_size/2
binding_sites_U21$start = binding_sites_U21$position_effect_size - binding_sites_U21$width_effect_size/2
binding_sites_U21$end = binding_sites_U21$position_effect_size + binding_sites_U21$width_effect_size/2

### TF -----
binding_sites_C03_ranges = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_L12_ranges = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_Q17_ranges = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)
binding_sites_U21_ranges = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width_effect_size <= 70, ], 
                                                    keep.extra.columns = T)

overlaps_loci_c03 = findOverlaps(aging_loci, binding_sites_C03_ranges)
overlaps_loci_l12 = findOverlaps(aging_loci, binding_sites_L12_ranges)
overlaps_loci_q17 = findOverlaps(aging_loci, binding_sites_Q17_ranges)
overlaps_loci_u21 = findOverlaps(aging_loci, binding_sites_U21_ranges)

loci_c03 = aging_loci[overlaps_loci_c03@from]
loci_l12 = aging_loci[overlaps_loci_l12@from]
loci_q17 = aging_loci[overlaps_loci_q17@from]
loci_u21 = aging_loci[overlaps_loci_u21@from]
c(length(unique(loci_c03)), 
  length(unique(loci_l12)), 
  length(unique(loci_q17)), 
  length(unique(loci_u21)))

### nucleosome -----
binding_sites_C03_nucle = makeGRangesFromDataFrame(binding_sites_C03[binding_sites_C03$width > 120, ])
binding_sites_L12_nucle = makeGRangesFromDataFrame(binding_sites_L12[binding_sites_L12$width > 120, ])
binding_sites_Q17_nucle = makeGRangesFromDataFrame(binding_sites_Q17[binding_sites_Q17$width > 120, ])
binding_sites_U21_nucle = makeGRangesFromDataFrame(binding_sites_U21[binding_sites_U21$width > 120, ])

overlaps_loci_c03_nucle = findOverlaps(aging_loci, binding_sites_C03_nucle)
overlaps_loci_l12_nucle = findOverlaps(aging_loci, binding_sites_L12_nucle)
overlaps_loci_q17_nucle = findOverlaps(aging_loci, binding_sites_Q17_nucle)
overlaps_loci_u21_nucle = findOverlaps(aging_loci, binding_sites_U21_nucle)
c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
  length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from))

loci_c03_nucle = aging_loci[overlaps_loci_c03_nucle@from]
loci_l12_nucle = aging_loci[overlaps_loci_l12_nucle@from]
loci_q17_nucle = aging_loci[overlaps_loci_q17_nucle@from]
loci_u21_nucle = aging_loci[overlaps_loci_u21_nucle@from]

length(findOverlaps(loci_c03_nucle, loci_l12_nucle)@from)/length(loci_c03_nucle)
length(findOverlaps(loci_l12_nucle, loci_q17_nucle)@from)/length(loci_l12_nucle)
length(findOverlaps(loci_q17_nucle, loci_u21_nucle)@from)/length(loci_q17_nucle)

overlaps_loci = data.frame('age' = c("weeks_16", 
                                     "weeks_30", 
                                     "weeks_56", 
                                     "weeks_82"), 
                           'TF' = c(length(overlaps_loci_c03@from), length(overlaps_loci_l12@from), 
                                    length(overlaps_loci_q17@from), length(overlaps_loci_u21@from)), 
                           'Nucleosome' = c(length(overlaps_loci_c03_nucle@from), length(overlaps_loci_l12_nucle@from), 
                                            length(overlaps_loci_q17_nucle@from), length(overlaps_loci_u21_nucle@from)))
overlaps_loci = overlaps_loci %>%
  pivot_longer(cols = c(TF, Nucleosome), 
               values_to = 'count', 
               names_to = 'Footprinting')
p1 = ggplot(overlaps_loci) + 
  geom_bar(aes(x = Footprinting, y = count, fill = age), 
           stat="identity", 
           position="dodge") + 
  ylab("counts of overlapped age-associated loci") + 
  scale_fill_manual(values = c('#D3E671',
                               '#A3D8FF', 
                               '#FAA301',
                               '#E58C94'), 
                    labels = c("weeks_16" = "16 weeks", 
                               "weeks_30" = "30 weeks", 
                               "weeks_56" = "56 weeks", 
                               "weeks_82" = "82 weeks")) +
  theme_Publication()
p1
ggsave("counts_overlapped_age-associated_loci.png", 
       p1, 
       height = 5, width = 7)

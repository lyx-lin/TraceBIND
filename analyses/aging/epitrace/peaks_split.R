library(tidyverse)
library(GenomicRanges)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all"))
# source("../../../code/utils.R")
# source("../../../code/getBias.R")
source("~/codes/footprints/footprint_prediction.R")

regionsBed = read.table("all_peaks_region.bed")
colnames(regionsBed) = c('chr', 'start', 'end')

### binding sites ====
counts_16 = readRDS("age_16_weeks/C03/counts.rds")
counts_30 = readRDS("age_30_weeks/L12/counts.rds")
counts_56 = readRDS("age_56_weeks/Q17/counts.rds")
counts_82 = readRDS("age_82_weeks/U21/counts.rds")

logp_16 = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_16_weeks/C03/kidney_multiome/logp_mean_threshold.txt"))
logp_30 = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_30_weeks/L12/kidney_multiome/logp_mean_threshold.txt"))
logp_56 = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_56_weeks/Q17/kidney_multiome/logp_mean_threshold.txt"))
logp_82 = read.table(paste0("~/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND/age_82_weeks/U21/kidney_multiome/logp_mean_threshold.txt"))

logp = logp_16
logp$smoothed_threshold = apply(cbind(logp_16$smoothed_threshold, 
                                      logp_30$smoothed_threshold, 
                                      logp_56$smoothed_threshold, 
                                      logp_82$smoothed_threshold), 1, mean)
logp$threshold = apply(cbind(logp_16$threshold, 
                             logp_30$threshold, 
                             logp_56$threshold, 
                             logp_82$threshold), 1, mean)

setwd("~/nzhanglab/project/linyx/footprints/results/data/Parker_DEFND_aging/all/all_ages/effect_size")
files = list.files()
files = gtools::mixedsort(files)
binding_sites_all_thresholding = pbmcapply::pbmclapply(1:length(files), 
                                                       function(i) {
                                                         binding_sites = read.table(files[i], header = TRUE) 
                                                         j = as.numeric(strsplit(files[i], '_')[[1]][1])
                                                         if (dim(binding_sites)[1] > 0){
                                                           binding_sites$threshold = logp$smoothed_threshold[which.min(abs(as.numeric(rownames(logp)) - binding_sites$mean_coverage))]
                                                           binding_sites$chr <- rep(regionsBed$chr[j], dim(binding_sites)[1])
                                                           binding_sites = binding_sites[binding_sites$p_value >= binding_sites$threshold, ]
                                                         }
                                                         binding_sites
                                                       }, 
                                                       mc.cores = 2
)
saveRDS(binding_sites_all_thresholding, 
        '/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all/thresholded_binding_sites_PT.rds')

### merging seuratobject ----
library(Signac)
library(Seurat)
metadata = read.table("~/nzhanglab/data/ParkerWilson_DEFND/metadata.txt")
cells_conversion = read.csv("~/nzhanglab/data/ParkerWilson_DEFND/barcode_conversion_10Xmultiome.csv")
### 16 weeks ----
metadata_16 = metadata[metadata$library_id == "RAGE24-C03-KYC-LN-01-Multiome-RNA", ]
metadata_16$GEX_bc = rownames(metadata_16)
metadata_16$GEX_bc = gsub(substr(metadata_16$GEX_bc[1], nchar(metadata_16$GEX_bc[1]) + 1 - 2, nchar(metadata_16$GEX_bc[1])), '', metadata_16$GEX_bc)

barcodeGroups_16 = left_join(metadata_16, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_16$GEX_bc = NULL
barcodeGroups_16$ATAC_bc = paste0(barcodeGroups_16$ATAC_bc, '-1')
frags.16 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-C03-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_16$ATAC_bc
)

### 30 weeks -----
metadata_30 = metadata[metadata$library_id == "RAGE24-L12-KYC-LN-01-Multiome-RNA", ]
metadata_30$GEX_bc = rownames(metadata_30)
metadata_30$GEX_bc = gsub(substr(metadata_30$GEX_bc[1], nchar(metadata_30$GEX_bc[1]) + 1 - 3, nchar(metadata_30$GEX_bc[1])), '', metadata_30$GEX_bc)

barcodeGroups_30 = left_join(metadata_30, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_30$GEX_bc = NULL
barcodeGroups_30$ATAC_bc = paste0(barcodeGroups_30$ATAC_bc, '-1')
frags.30 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-L12-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_30$ATAC_bc
)

### 56 weeks -----
metadata_56 = metadata[metadata$library_id == "RAGE24-Q17-KYC-LN-01-Multiome-RNA", ]
metadata_56$GEX_bc = rownames(metadata_56)
metadata_56$GEX_bc = gsub(substr(metadata_56$GEX_bc[1], nchar(metadata_56$GEX_bc[1]) + 1 - 3, nchar(metadata_56$GEX_bc[1])), '', metadata_56$GEX_bc)

barcodeGroups_56 = left_join(metadata_56, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_56$GEX_bc = NULL
barcodeGroups_56$ATAC_bc = paste0(barcodeGroups_56$ATAC_bc, '-1')
frags.56 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-Q17-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_56$ATAC_bc
)

### 82 weeks ----
metadata_82 = metadata[metadata$library_id == "RAGE24-U21-KYC-LN-01-Multiome-RNA", ]
metadata_82$GEX_bc = rownames(metadata_82)
metadata_82$GEX_bc = gsub(substr(metadata_82$GEX_bc[1], nchar(metadata_82$GEX_bc[1]) + 1 - 3, nchar(metadata_82$GEX_bc[1])), '', metadata_82$GEX_bc)

barcodeGroups_82 = left_join(metadata_82, 
                             cells_conversion, by = 'GEX_bc')
barcodeGroups_82$GEX_bc = NULL
barcodeGroups_82$ATAC_bc = paste0(barcodeGroups_82$ATAC_bc, '-1')
frags.82 <- CreateFragmentObject(
  path = "~/nzhanglab/data/ParkerWilson_DEFND/RAGE24-U21-KYC-LN-01-Multiome-ATAC/fragments.tsv.gz",
  cells = barcodeGroups_82$ATAC_bc
)

### split peaks ====
setwd("/home/mnt/weka/nzh/nzhanglab/project/linyx/footprints/PRINT/data/Parker_DEFND_aging/all")
regionsBed = read.table("all_peaks_region.bed")
peaks = makeGRangesFromDataFrame(regionsBed)

binding_sites_all = readRDS("thresholded_binding_sites_PT.rds")
binding_sites_all = data.table::rbindlist(binding_sites_all, fill = T)
binding_sites_all$start = binding_sites_all$position - pmax(50, binding_sites_all$width*3/2)
binding_sites_all$end = binding_sites_all$position + pmax(50, binding_sites_all$width*3/2)
split_gr = makeGRangesFromDataFrame(binding_sites_all, 
                                                    keep.extra.columns = T)
split_gr = reduce(split_gr)
split_gr <- disjoin(c(peaks, split_gr))

### combining -----
weeks16.counts <- FeatureMatrix(
  fragments = frags.16,
  features = split_gr,
  cells = barcodeGroups_16$ATAC_bc
)

weeks30.counts <- FeatureMatrix(
  fragments = frags.30,
  features = split_gr,
  cells = barcodeGroups_30$ATAC_bc
)

weeks56.counts <- FeatureMatrix(
  fragments = frags.56,
  features = split_gr,
  cells = barcodeGroups_56$ATAC_bc
)

weeks82.counts <- FeatureMatrix(
  fragments = frags.82,
  features = split_gr,
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

seuratobject_fp <- merge(
  x = weeks16,
  y = list(weeks30, weeks56, weeks82),
  add.cell.ids = c("weeks_16", "weeks_30", "weeks_56", "weeks_82")
)

saveRDS(seuratobject_fp, 
        'seuratobject_merged_split_peaks_thresholded_allcelltypes.rds')

### visualization ----
library(Signac)
library(Seurat)
seuratobject_fp = readRDS("seuratobject_merged_split_peaks_no_smoothed_PT_thresholded_95_effect_size.rds")
seuratobject_fp <- RunTFIDF(seuratobject_fp)
seuratobject_fp <- FindTopFeatures(seuratobject_fp, 
                                      min.cutoff = 20)
seuratobject_fp <- RunSVD(seuratobject_fp)

seuratobject_fp <- harmony::RunHarmony(
  object = seuratobject_fp,
  group.by.vars = 'dataset',
  assay.use = 'ATAC',
  reduction.use = 'lsi',
  reduction.save = 'harmony',
  project.dim = FALSE
)

seuratobject_fp <- RunUMAP(seuratobject_fp, 
                              dims = 2:30, 
                              reduction = 'harmony', 
                              reduction.name = 'umap_harmony')

beach = c("pink", '#99CC66', "#00A08A", "#0B775E", 
          "#87D2DB","#5BB1CB",'#88CEEF', '#1E78B4', "#4F66AF", "#352A86", 
          '#a8ddb5',"#FF6800","#FCBF6E", "#fbdf72", 
          "#8c6bb1", "#7F3F98", "#FCB31A", '#E6C2DC')
p2 = ggplot() + 
  geom_point(aes(x = Embeddings(seuratobject_fp[['umap_harmony']])[, 1], 
                 y = Embeddings(seuratobject_fp[['umap_harmony']])[, 2], 
                 color = seuratobject_fp$celltype), 
             size = 1) +  
  # ggrepel::geom_text_repel(data = centroids_harmony_split, 
  #                          aes(x = umap_1, y = umap_2, label = celltype),
  #                          size = 4, color = "black") +  # Add labels
  scale_colour_manual(values = beach, 
                      name = 'cell type') + 
  ylab("umap_2") + 
  xlab("umap_1") + 
  theme_Publication()
p2
ggsave("umap_integration_clustering_splitup_peaks_celltypes.pdf", 
       height = 8, width = 7)
summerNight = c("#FCBF6E", "#bed678","#CAB2D6","#A1CDE1")
p2 = ggplot() + 
  geom_point(aes(x = Embeddings(seuratobject_fp[['umap_harmony']])[, 1], 
                 y = Embeddings(seuratobject_fp[['umap_harmony']])[, 2], 
                 color = seuratobject_fp$dataset), 
             size = 1) + 
  scale_colour_manual(values = summerNight, 
                      name = 'ages') + 
  ylab("umap_2") + 
  xlab("umap_1") + 
  theme_Publication()
ggsave("umap_integration_clustering_splitup_peaks_samples.pdf", 
       height = 6, width = 7)


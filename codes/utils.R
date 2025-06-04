library(reticulate)
library(MASS)
library(tidyverse)
library(GenomicRanges)

### get insertions from the fragments from the group groupIDs with width
get_insertions <- function(countData, # Tn5 insertion count tensor
                           regionInd, # Index of the region
                           groupIDs, 
                           width # Width of the region
){
  if (dim(countData[[regionInd]])[1] == 0){
    regionTrack <- rep(0, width)
    regionTrack
  } else {
    regionTracks = sapply(
      groupIDs,
      function(groupID){
        groupRegionATAC <- countData[[regionInd]] %>% filter(group %in% groupID)
        regionTrack <- rep(0, width)
        regionTrack[groupRegionATAC$position] <- groupRegionATAC$count
        regionTrack
      }
    )
    apply(data.frame(regionTracks), 1, sum)
  }
}

get_count = function(frags_path, 
                     regionsBed, 
                     barcodeGroups = NULL, 
                     chunkSize = 2000, 
                     tmpDir = "./", 
                     nCores = 2, 
                     pairedEnd = TRUE, 
                     startsAre0based = TRUE, 
                     bulk = FALSE){
  frags = data.table::fread(frags_path, 
                            # sep = '\t', 
                            showProgress = TRUE) %>% 
    as.data.frame()
  
  colnames(frags)[1:3] = c("V1", "V2", "V3")
  
  # colnames(frags) = c("V1", "V2", "V3", "V4", 'V5')
  
  # frags = apply(regionsBed, 1, function(x) {
  #   frags[(frags$V1 == x[1] & 
  #            frags$V2 <= as.numeric(x[2]) & 
  #            frags$V3 >= as.numeric(x[2]))|
  #           (frags$V1 == x[1] & 
  #              frags$V2 <= as.numeric(x[3]) & 
  #              frags$V3 >= as.numeric(x[3]))| 
  #           (frags$V1 == x[1] & 
  #              frags$V2 >= as.numeric(x[2]) & 
  #              frags$V3 <= as.numeric(x[3])), ]
  # })
  
  # frags = data.table::rbindlist(frags)
  
  regions = GenomicRanges::GRanges(seqnames = regionsBed$chr, 
                                   ranges = IRanges::IRanges(start = regionsBed$start, 
                                                             end = regionsBed$end))
  # width(regions) = max(width(regions)) 
  if (is.null(barcodeGroups)){
    barcodeGroups = unique(as.data.frame(frags)[4])
    barcodeGroups['group'] = 1
  } 
  
  if (dim(barcodeGroups)[1] == 1){
    barcodeGroups['group'] = 1
  }
  colnames(barcodeGroups) = c('barcodeID', 'group')
  
  frags = frags %>% GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", 
                                                            start.field = "V2", 
                                                            end.field = "V3", 
                                                            keep.extra.columns = TRUE, 
                                                            starts.in.df.are.0based = startsAre0based)
  
  counts = computeCountTensor(frags, 
                              regions, 
                              barcodeGroups, 
                              chunkSize = chunkSize, 
                              pairedEnd = pairedEnd, 
                              nCores = nCores, 
                              bulk = bulk)
}

computeCountTensor <- function(frags, # Path to fragments file
                               regions, # Genomic ranges of the regions to footprint
                               barcodeGroups, # Data.frame specifying membership of barcodes in pseudobulks. First column is barcodes and second is groupID 
                               chunkSize = 2000, # Chunk size for parallel processing of regions
                               nCores = 2, # Number of cores to use
                               pairedEnd = TRUE, # # True/False. Whether data is paired-end or single-end
                               bulk = FALSE
                               # returnCombined = TRUE, # Whether to return the combined result for all chunks. Set it to False when data is too big 
                               # file_path = NULL
) {
  start_time <- Sys.time()
  
  #################### Load fragments data, convert to GRanges #######################
  
  cat("Make 1bp step .. \n")
  nRegions <- length(regions)
  
  # Rename extra columns
  if (ncol(mcols(frags)) > 1) {
    colnames(mcols(frags)) <- c("barcodeID", "pcrDup")
  } else {
    colnames(mcols(frags)) <- "barcodeID"
  }
  
  # Filter fragments by size (optional)
  
  # Get group ID for each cell
  groupInds <- gtools::mixedsort(unique(barcodeGroups$group))
  gc()
  
  #################### Get region-by-position-by-pseudobulk Tn5 insertion count tensor #######################
  
  # To reduce memory usage, we chunk the data in to smaller chunks
  cat("Reformating counts data into a list (each element is data for a region) ..\n")
  chunkIntervals <- getChunkInterval(regions, chunkSize = chunkSize)
  starts <- chunkIntervals[["starts"]]
  ends <- chunkIntervals[["ends"]]
  
  # Create a folder for saving intermediate results
  
  # For each chunk, we extract data for individual regions
  # For each region, we store the data in a 3-column data.frame (columns are group, position and counts)
  # Re-organize data into lists
  countTensorAll = List()
  for(i in 1:length(starts)){
    
    print("Re-organizing data into lists")
    print(paste0(Sys.time()," Processing chunk ", i, " out of ", length(starts), " chunks"))
    
    # Skip current chunk if result already exists
    
    # Get fragments within the current chunk
    chunkRegions <- starts[i]:ends[i]
    chunkFrags <- subsetByOverlaps(frags, regions[chunkRegions])
    
    # Go through each group (pseudobulk) and retrive cutsite data
    # Should result in a data table with 4 columns: Region index, position in this region, group ID, count
    groupedCountTensor <- pbmcapply::pbmclapply(
      groupInds,
      function(groupInd){
        
        # Get fragmenst belonging to the current cell group (i.e. pseudobulk)
        groupBarcodes <- barcodeGroups$barcodeID[barcodeGroups$group %in% groupInd]
        groupFrags <- chunkFrags[chunkFrags$barcodeID %in% groupBarcodes]
        
        if(pairedEnd){
          # Get all Tn5 insertion sites (single base pair resolution)
          # Note: The input fragments file should be +4/-5 shifted to accommodate common practice
          # However, +4/-5 actually points to the base immediately to the left of the center of the 9bp staggered end
          # The 1bp cut position should actually by +5/-4. Therefore we need to shift both start and end by +1
          # The function fragsToRanges already shifts start position by +1 when we specify "startsAre0based = T"
          # Therefore here we only need to further shift the end position by + 1
          cutsites <- c(resize(groupFrags, width = 1, fix = "start"),
                        IRanges::shift(resize(groupFrags, width = 1, fix = "end"), 1))
        } else {
          # If we were provided single-end data, the fragments file should be in the format of
          # (chr name, cut position, cut position, barcode, number of insertions at this pos)
          cutsites <- c(resize(groupFrags, width = 1, fix = "start"))
        }
        
        # Get position of cutsites within regions
        ovRegions <- findOverlaps(query = regions, 
                                  subject = cutsites)
        if(length(ovRegions) == 0){
          groupCountTensor <- NULL
        } else {
          positions <- start(cutsites)[ovRegions@to] - start(regions)[ovRegions@from] + 1
          # Generate a data frame with 4 columns: Region index, position in this region, group ID, count
          if (bulk){
            x = cutsites[ovRegions@to]$pcrDup
          } else {
            x = 1
          }
          cat(paste0("Generating matrix of counts for group ", groupInd, "..\n"))
          groupCountTensor <- Matrix::sparseMatrix(i = ovRegions@from,
                                                   j = positions,
                                                   x = x)
          groupCountTensor <- summary(t(groupCountTensor))
          colnames(groupCountTensor) <- c("position", "region", "count")
          groupCountTensor$group <- groupInd
          groupCountTensor <- groupCountTensor[,c("region", "position", "group", "count")]
          groupCountTensor <- as_tibble(groupCountTensor)
        }
        
        groupCountTensor
      }, 
      mc.cores = nCores
    )
    groupedCountTensor <- data.table::rbindlist(groupedCountTensor)
    
    # Re-organize the data into a list where each element is the data for a single region.
    if(dim(groupedCountTensor)[1] > 0){
      cluster <- prep_cluster(length(chunkRegions), n_cores = nCores)
      opts <- cluster[["opts"]]
      cl <- cluster[["cl"]]
      countTensorChunk <- foreach(regionInd = chunkRegions,
                                  .options.snow = opts, 
                                  .packages = c("dplyr","Matrix")) %dopar%   {
                                    regionCountTensor <- groupedCountTensor %>% filter(region %in% regionInd)
                                    return(regionCountTensor)
                                  }
      stopCluster(cl)
    }else{
      countTensorChunk <- lapply(chunkRegions, function(x){data.table::data.table()})
    }
    
    countTensorAll = c(countTensorAll, countTensorChunk)
    
    # Save results
    # Release unused memory
    if((i %% 10) == 0) gc()
  }
  
  cat("Done!\n")
  end_time <- Sys.time()
  cat("Time elapsed: ", end_time - start_time, units(end_time - start_time), " \n\n")
  return(countTensorAll)
}

get_bias = function(regionsBed,
                    referenceGenome = 'hg38', 
                    genome = NULL, 
                    contextRadius = 50,
                    nCores = 2,
                    chunkSize = 2000, 
                    model_use = 'Tn5_NN_model.h5', 
                    path = NULL, 
                    code_path = '../codes/', 
                    model_path = '../../'){
  if (is.null(path)){
    path = getwd()
  }
  
  if(referenceGenome == "hg19"){
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }else if(referenceGenome == "hg38"){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else if(referenceGenome == "mm10"){
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }else if(referenceGenome == "rn7"){
    genome <- BSgenome.Rnorvegicus.UCSC.rn7::BSgenome.Rnorvegicus.UCSC.rn7
  }
  
  regions = GenomicRanges::GRanges(seqnames = regionsBed$chr, 
                                   ranges = IRanges::IRanges(start = regionsBed$start, 
                                                             end = regionsBed$end))
  width(regions) = max(width(regions))
  contextLen = 2*contextRadius + 1
  if(!file.exists(paste(path, "regionSeqs.txt", sep = ""))){
    regionSeqs = pbmcapply::pbmclapply(
      1:length(regions),
      function(regionInd){
        range = regions[regionInd]
        # Get genomic sequence of the genomic region
        # We extract an extra flanking region of length contextLen on both sides (so we can predict bias for edge positions)
        regionSeq = Biostrings::getSeq(genome, as.character(seqnames(range)), 
                                       start = start(range) - contextRadius, 
                                       width = width(range) + contextLen - 1,
                                       as.character = T)
        regionSeq
      },
      mc.cores = nCores
    )
    # Save sequence context to a file
    regionSeqs = as.character(regionSeqs)
    write.table(regionSeqs, paste(path, "regionSeqs.txt", sep = ""), quote = F,
                col.names = F, row.names = F)
  }
  
  write.table(c(model_path, path, model_use, chunkSize), 
              "args.txt", quote = F, col.names = F, row.names = F)
  
  py_run_file(paste0(code_path, "/predictBias.py"))
  # py_run_file(paste0("/home/mnt/weka/nzh/team/woodsqu2/nzhanglab/project/linyx/footprints/PRINT/code/predictBias.py"))

  pred_bias = read.table(paste0(path, "pred_bias.txt"))
  pred_bias[pred_bias <= 0] = 1e-10
  write.table(pred_bias, paste0(path, "pred_bias.txt"))
}

finetuned_model = function(code_path, 
                           obsbias_path, 
                           PRINT_model_path, 
                           finetuned_model_save_path, 
                           finetuned_model_name
                           ){
  write.table(c(obsbias_path, 
                PRINT_model_path, 
                finetuned_model_save_path, 
                finetuned_model_name), 
              "args.txt", quote = F, col.names = F, row.names = F)
  py_run_file(paste0(code_path, "/finetuning.py"))
}
  
## create the bias and insertions data for the peak
get_bardata = function(counts, 
                       regionsBed, 
                       bias, 
                       i, 
                       groupIDs){
  insertions = get_insertions(counts, 
                              i, 
                              groupIDs = groupIDs, 
                              width = length(regionsBed$start[i]:regionsBed$end[i]))
  pred_bias = as.numeric(bias)
  bardata = data.frame("position" = as.numeric(format(regionsBed$start[i]:regionsBed$end[i], scientific = FALSE, trim = TRUE)),
                       "Tn5Insertion" = insertions,
                       "pred_bias" = pred_bias)
  bardata
}


mito_fdr = function(mito_counts, 
                    mito_bias, 
                    mito_regions, 
                    i, 
                    seeds, 
                    average_insertions, 
                    alpha = 0.05, 
                    ncore = 2){
  thresholds = lapply(seeds, function(k) {
    set.seed(k)
    mito_downsamples(mito_counts, 
                     mito_bias, 
                     mito_regions, 
                     i, 
                     average_insertions, 
                     ncore = ncore)
    }
  )
  thresholds = do.call(rbind, thresholds)
  labels = seq(0, 40, 0.2)
  breaks <- c(0, 
              sapply(1:(length(labels) - 1), function(i) mean(c(labels[i], labels[i + 1]))), Inf)
  labels <- as.character(labels)  # Corresponding labels
  thresholds$labels <- cut(as.numeric(thresholds$mean_coverage), 
                                       breaks = breaks, 
                                       labels = labels, 
                                       right = F)
  logp = thresholds %>%
    group_by(labels) %>%
    summarise(
      count = n(), 
      # label1 = as.numeric(unique(labels)), 
      threshold = quantile(p_value, alpha)
      # threshold = mean(p_value) + 2.32*sd(p_value)/n()
      # threshold = median(p_value) + 1.57*IQR(p_value)/sqrt(1 + n())
    )
  logp = na.omit(logp)
  logp$smoothed_threshold = logp$threshold
  logp$smoothed_threshold[3:(dim(logp)[1] - 2)] = sapply(3:(dim(logp)[1]-2), 
                                                         function(i) mean(logp$threshold[(i-2):(i+2)]))
  logp = logp[1:150, ]
  return(logp)
}

mito_downsample = function(mito_counts, 
                           mito_bias, 
                           mito_regions, 
                           i, 
                           average_insertion, 
                           ncore = 2){
  bardata = get_bardata(counts = mito_counts, 
                        regionsBed = mito_regions, 
                        bias = mito_bias[i, 1:length(mito_regions$start[i]:mito_regions$end[i])], 
                        i = i, 
                        groupIDs = unique(mito_counts[[i]]$group))
  ratio = average_insertion/mean(bardata$Tn5Insertion)
  bardata$downsampled_Tn5Insertion = apply(bardata['Tn5Insertion'], 1, function(x) 
    rbinom(n = 1, size = x, prob = ratio))
  
  footprinting_results = NB_footprintings(Tn5Insertion = bardata$downsampled_Tn5Insertion, 
                                          pred_bias = bardata$pred_bias, 
                                          positions = bardata$position, 
                                          p.adjust.method = 'BH', 
                                          nCores = ncore)
  p_value_matrix = footprinting_results[['pval']]
  effect_size_matrix = footprinting_results[['effect_size']]
  
  binding_sites_freezed = binding_sites_pval(p_value_matrix, 
                                             width_threshold = 10, 
                                             pval_threshold = 1)
  
  binding_sites_freezed['mean_coverage'] = rep(NA, dim(binding_sites_freezed)[1])
  
  if (dim(binding_sites_freezed)[1] > 0){
    for (index in (1:dim(binding_sites_freezed)[1])){
      positions_1 = (binding_sites_freezed$position[index] - binding_sites_freezed$width[index]/2 - 1):(binding_sites_freezed$position[index] - 1 - binding_sites_freezed$width[index]*3/2)
      positions_2 = (binding_sites_freezed$position[index] + binding_sites_freezed$width[index]/2 + 1):(binding_sites_freezed$position[index] + 1 + binding_sites_freezed$width[index]*3/2)
      mean_coverage = min(mean(bardata[bardata$position %in% positions_1, 'downsampled_Tn5Insertion']), 
                          mean(bardata[bardata$position %in% positions_2, 'downsampled_Tn5Insertion']))
      
      binding_sites_freezed[index, 'mean_coverage'] = mean_coverage
    }
  }
  return(binding_sites_freezed)
}

mito_downsamples = function(mito_counts, 
                            mito_bias, 
                            mito_regions, 
                            i, 
                            average_insertions, 
                            ncore = 2){
  binding_sites = lapply(average_insertions, 
                         function(average_insertion){ 
                           mito_downsample(mito_counts, 
                                           mito_bias, 
                                           mito_regions, 
                                           i=i, 
                                           average_insertion=average_insertion, 
                                           ncore = ncore)
                         }
  )
  binding_sites <- do.call(rbind, binding_sites)
  return(binding_sites)
}
getChunkInterval <- function(x, # Vector or list
                             chunkSize = 2000 # Size of a single chunk
){
  
  chunkSize <- min(length(x), chunkSize)
  nData <- length(x)
  starts <- seq(1, nData, chunkSize)
  ends <- starts + chunkSize - 1
  ends[length(ends)] <- nData
  
  list("starts" = starts,
       "ends" = ends)
}

prep_cluster <- function(len, # Number of elemnts in the iterable list
                         n_cores = 2 # Number of cores to use
){
  library(doParallel)
  opts <- list()
  pb <- txtProgressBar(min = 0, max = len, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- parallel::makeCluster(n_cores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  list("opts" = opts, "cl" = cl)
}
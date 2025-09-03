#' @title get_insertions
#'
#' @description For given group, get insertions
#' @details For given group, get insertions from count data on a particular region
#' @param countData Tn5 insertion count on a particular region from get_count
#' @param groupIDs Group name
#' @param width Width of the region
#' @import dplyr
#' @return A vector containing insertions for every base pair
#' @export
get_insertions <- function(countData,
                           regionInd,
                           groupIDs,
                           width
){
  if (dim(countData[[regionInd]])[1] == 0){
    regionTrack <- rep(0, width)
    regionTrack
  } else {
    regionTracks = sapply(
      groupIDs,
      function(groupID){
        groupRegionATAC <- countData[[regionInd]] %>% dplyr::filter(group %in% groupID)
        regionTrack <- rep(0, width)
        regionTrack[groupRegionATAC$position] <- groupRegionATAC$count
        regionTrack
      }
    )
    apply(data.frame(regionTracks), 1, sum)
  }
}

#' @title get_count
#'
#' @description get insertion count for each base pair
#' @details get insertion count for each base pair from frag file
#' @param frags_path path to frags
#' @param regionsBed a data.frame for peak, containing columns chr, start. end
#' @param barcodeGroups a data.frame for barcodes in pseudobulks. First column is barcodes and second is group info
#' @param chunkSize Chunk size for parallel processing of regions
#' @param nCores Number of cores to use
#' @param pairedEnd Whether data is paired-end or single-end
#' @param startsAre0based Whether the coordinates are 0-based
#' @import GenomicRanges
#' @import IRanges
#' @importFrom data.table fread
#' @return A count tensor file
#' @export
get_count = function(frags_path,
                     regionsBed,
                     barcodeGroups = NULL,
                     chunkSize = 2000,
                     nCores = 2,
                     pairedEnd = TRUE,
                     startsAre0based = TRUE){
  frags = fread(frags_path,
                showProgress = TRUE) %>%
    as.data.frame()

  colnames(frags)[1:3] = c("V1", "V2", "V3")

  regions = GRanges(seqnames = regionsBed$chr,
                    ranges = IRanges(start = regionsBed$start,
                                     end = regionsBed$end))

  if (is.null(barcodeGroups)){
    barcodeGroups = unique(as.data.frame(frags)[4])
    barcodeGroups['group'] = 1
  }

  if (dim(barcodeGroups)[1] == 1){
    barcodeGroups['group'] = 1
  }
  colnames(barcodeGroups) = c('barcodeID', 'group')

  frags = frags %>% makeGRangesFromDataFrame(seqnames.field = "V1",
                                             start.field = "V2",
                                             end.field = "V3",
                                             keep.extra.columns = TRUE,
                                             starts.in.df.are.0based = startsAre0based)

  counts = computeCountTensor(frags,
                              regions,
                              barcodeGroups,
                              chunkSize = chunkSize,
                              pairedEnd = pairedEnd,
                              nCores = nCores)
}

#' @title computeCountTensor
#'
#' @description help function for get insertion count for each base pair
#' @details help function for get insertion count for each base pair from frag file
#' @param frags_path path to frags
#' @param regions Genomic ranges of the regions to footprint
#' @param barcodeGroups Data.frame specifying membership of barcodes in pseudobulks. First column is barcodes and second is groupID
#' @param chunkSize Chunk size for parallel processing of regions
#' @param nCores Number of cores to use
#' @param pairedEnd Whether data is paired-end or single-end
#' @import dplyr
#' @import GenomicRanges
#' @import IRanges
#' @import Matrix
#' @import pbmcapply
#' @import tibble
#' @importFrom data.table rbindlist
#' @importFrom foreach foreach
#' @importFrom gtools mixedsort
#' @return A count tensor file
computeCountTensor <- function(frags,
                               regions,
                               barcodeGroups,
                               chunkSize = 2000,
                               nCores = 2,
                               pairedEnd = TRUE
) {
  start_time <- Sys.time()

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
  groupInds <- mixedsort(unique(barcodeGroups$group))
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
    groupedCountTensor <- pbmclapply(
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
          cutsites <- c(GenomicRanges::resize(groupFrags, width = 1, fix = "start"),
                        GenomicRanges::shift(GenomicRanges::resize(groupFrags, width = 1, fix = "end"), 1))
        } else {
          # If we were provided single-end data, the fragments file should be in the format of
          # (chr name, cut position, cut position, barcode, number of insertions at this pos)
          cutsites <- c(GenomicRanges::resize(groupFrags, width = 1, fix = "start"))
        }

        # Get position of cutsites within regions
        ovRegions <- findOverlaps(query = regions,
                                  subject = cutsites)
        if(length(ovRegions) == 0){
          groupCountTensor <- NULL
        } else {
          positions <- GenomicRanges::start(cutsites)[ovRegions@to] - GenomicRanges::start(regions)[ovRegions@from] + 1
          # Generate a data frame with 4 columns: Region index, position in this region, group ID, count
          cat(paste0("Generating matrix of counts for group ", groupInd, "..\n"))
          groupCountTensor <- Matrix::sparseMatrix(i = ovRegions@from,
                                                   j = positions,
                                                   x = 1)
          groupCountTensor <- Matrix::summary(Matrix::t(groupCountTensor))
          colnames(groupCountTensor) <- c("position", "region", "count")
          groupCountTensor$group <- groupInd
          groupCountTensor <- groupCountTensor[, c("region", "position", "group", "count")]
          groupCountTensor = as_tibble(groupCountTensor)
        }
        groupCountTensor
      },
      mc.cores = nCores
    )
    groupedCountTensor <- rbindlist(groupedCountTensor)

    if(dim(groupedCountTensor)[1] > 0){
      cluster <- prep_cluster(length(chunkRegions), n_cores = nCores)
      opts <- cluster[["opts"]]
      cl <- cluster[["cl"]]
      countTensorChunk <- foreach(regionInd = chunkRegions,
                                  .options.snow = opts,
                                  .packages = c("dplyr","Matrix")) %dopar%   {
                                    regionCountTensor <- groupedCountTensor %>% dplyr::filter(region %in% regionInd)
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

#' @title get_bias
#'
#' @description Get bias for each base pair
#' @details Get bias for each base pair by finetuned/trained model
#' @param regionsBed a data frame for peak, containing columns chr, start. end
#' @param referenceGenome hg19/hg38/mm10/rn7
#' @param genome Full genome sequences if not hg19/hg38/mm10/rn7
#' @param chunkSize Chunk size for parallel processing of regions
#' @param nCores Number of cores to use
#' @param model_use which bias model to use
#' @param model_path the path to model
#' @param save_path the path to save the bias file
#' @import GenomicRanges
#' @import IRanges
#' @import pbmcapply
#' @import reticulate
#' @importFrom BSgenome getSeq
#' @importFrom doSNOW registerDoSNOW
#' @return A bias tensor
#' @export
get_bias = function(regionsBed,
                    referenceGenome = 'hg38',
                    genome = NULL,
                    nCores = 2,
                    chunkSize = 2000,
                    model_use = 'Tn5_NN_model.h5',
                    model_path = '../../',
                    save_path = NULL){
  if (is.null(save_path)){
    save_path = getwd()
  }
  if(referenceGenome == "hg19"){
    if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
      BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }else if(referenceGenome == "hg38"){
    if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
      BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else if(referenceGenome == "mm10"){
    if (!require("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
      BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }else if(referenceGenome == "rn7"){
    if (!require("BSgenome.Rnorvegicus.UCSC.rn7", quietly = TRUE))
      BiocManager::install.packages("BSgenome.Rnorvegicus.UCSC.rn7")
    genome <- BSgenome.Rnorvegicus.UCSC.rn7::BSgenome.Rnorvegicus.UCSC.rn7
  } else {
    print('please specify your genome file in genome')
  }

  regions = GRanges(seqnames = regionsBed$chr,
                    ranges = IRanges(start = regionsBed$start,
                                     end = regionsBed$end))
  width(regions) = max(width(regions))
  contextRadius = 50
  contextLen = 2*contextRadius + 1
  if(!file.exists(paste(save_path, "regionSeqs.txt", sep = ""))){
    regionSeqs = pbmclapply(
      1:length(regions),
      function(regionInd){
        range = regions[regionInd]
        # Get genomic sequence of the genomic region
        # We extract an extra flanking region of length contextLen on both sides (so we can predict bias for edge positions)
        regionSeq = getSeq(genome, as.character(seqnames(range)),
                           start = start(range) - contextRadius,
                           width = width(range) + contextLen - 1,
                           as.character = T)
        regionSeq
      },
      mc.cores = nCores
    )
    # Save sequence context to a file
    regionSeqs = as.character(regionSeqs)
    write.table(regionSeqs, paste(save_path, "regionSeqs.txt", sep = ""), quote = F,
                col.names = F, row.names = F)
  }

  write.table(c(model_path, save_path, model_use, chunkSize),
              "args.txt", quote = F, col.names = F, row.names = F)

  py_run_file(system.file("python", "predictBias.py", package = "tracebind"))

  pred_bias = read.table(paste0(save_path, "pred_bias.txt"))
  pred_bias[pred_bias <= 0] = 1e-10
  write.table(pred_bias, paste0(save_path, "pred_bias.txt"))
}

#' @title finetuned_model
#'
#' @description Finetune PRINT model by observed bias
#' @details Finetune PRINT model by observed bias in mito DNA regions
#' @param obsbias_path path to file of observed bias in mito DNA regions
#' @param PRINT_model_path path to PRINT model
#' @param finetuned_model_save_path path to save the finetuned model
#' @param finetuned_model_name name to save the finetuned model
#' @param nCores Number of cores to use
#' @import reticulate
#' @return A finetuned model
#' @export
finetuned_model = function(obsbias_path,
                           PRINT_model_path,
                           finetuned_model_save_path,
                           finetuned_model_name,
                           nCores = 2
                           ){
  write.table(c(obsbias_path,
                PRINT_model_path,
                finetuned_model_save_path,
                finetuned_model_name,
                nCores),
              "args.txt", quote = F, col.names = F, row.names = F)
  py_run_file(system.file("python", "finetuning.py", package = "tracebind"))
}

#' @title get_bardata
#'
#' @description Create the bias and insertions data for the peak
#' @details create the bias and insertions data for the peak
#' @param counts the count tensor from get_count
#' @param regionsBed a data frame for peak, containing columns chr, start. end
#' @param bias the bias tensor from get_bias
#' @param i ID of peak of interest in regionsBed file
#' @param groupIDs group(s) of interest
#' @return A data frame containing position, insertions and predicted bias
#' @export
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

#' @title mito_fdr
#'
#' @description Create thresholds for False Discovery Control, stratified by average insertions
#' @details Create thresholds for stratified False Discovery Control based on mito regions
#' @param mito_counts the count tensor from get_count on mito regions
#' @param mito_bias the bias tensor from get_bias on mito regions
#' @param mito_regions a data frame for mito regions, containing columns chr, start. end
#' @param i ID of peak of interest in mito region, has to be the test region
#' @param seeds random seed for set.seed()
#' @param average_insertions average insertions of interest
#' @param alpha fdr
#' @param nCore Number of cores to use
#' @import dplyr
#' @return A data frame containing coverage, fdr threshold
#' @export
mito_fdr = function(mito_counts,
                    mito_bias,
                    mito_regions,
                    i,
                    seeds,
                    average_insertions,
                    alpha = 0.05,
                    nCore = 2){
  thresholds = lapply(seeds, function(k) {
    set.seed(k)
    mito_downsamples(mito_counts,
                     mito_bias,
                     mito_regions,
                     i,
                     average_insertions,
                     nCore = nCore)
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
    dplyr::group_by(labels) %>%
    dplyr::summarise(
      count = n(),
      threshold = quantile(p_value, alpha)
    )
  logp = na.omit(logp)
  logp$smoothed_threshold = logp$threshold
  logp$smoothed_threshold[3:(dim(logp)[1] - 2)] = sapply(3:(dim(logp)[1]-2),
                                                         function(i) mean(logp$threshold[(i-2):(i+2)]))
  logp = logp[1:150, ]
  return(logp)
}

#' @title mito_downsample
#'
#' @description Get footprints on downsampled mito DNA regions with a given average insertions
#' @details Get footprints on downsampled mito DNA regions with a given average insertions
#' @param mito_counts the count tensor from get_count on mito regions
#' @param mito_bias the bias tensor from get_bias on mito regions
#' @param mito_regions a data frame for mito regions, containing columns chr, start. end
#' @param i ID of peak of interest in mito region, has to be the test region
#' @param average_insertions average insertions of interest
#' @param nCore fdr
#' @return A data frame containing position, width, pvalue for non overlapping 'footprintings' (not necessarily significant)
#'
mito_downsample = function(mito_counts,
                           mito_bias,
                           mito_regions,
                           i,
                           average_insertion,
                           nCore = 2){
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

#' @title mito_downsamples
#'
#' @description Get footprints on downsampled mito DNA regions with all given average insertions
#' @details Get footprints on downsampled mito DNA regions with all given average insertions
#' @param mito_counts the count tensor from get_count on mito regions
#' @param mito_bias the bias tensor from get_bias on mito regions
#' @param mito_regions a data frame for mito regions, containing columns chr, start. end
#' @param i ID of peak of interest in mito region, has to be the test region
#' @param average_insertions average insertions of interest
#' @param nCore fdr
#' @return A data frame containing position, width, pvalue for non overlapping 'footprintings' (not necessarily significant)
#'
mito_downsamples = function(mito_counts,
                            mito_bias,
                            mito_regions,
                            i,
                            average_insertions,
                            nCore = 2){
  binding_sites = lapply(average_insertions,
                         function(average_insertion){
                           mito_downsample(mito_counts,
                                           mito_bias,
                                           mito_regions,
                                           i=i,
                                           average_insertion=average_insertion,
                                           nCore = nCore)
                         }
  )
  binding_sites <- do.call(rbind, binding_sites)
  return(binding_sites)
}

#' @title getChunkInterval
#'
#' @description Chunk a vector/list x into chunks. Return starts and ends of chunks
#' @details Chunk a vector/list x into chunks. Return starts and ends of chunks
#' @param x Vector or list
#' @param chunkSize Size of a single chunk
#' @import dplyr
#' @return start and end of x after chunking with chunkSize
#'
getChunkInterval <- function(x,
                             chunkSize = 2000
){

  chunkSize <- min(length(x), chunkSize)
  nData <- length(x)
  starts <- seq(1, nData, chunkSize)
  ends <- starts + chunkSize - 1
  ends[length(ends)] <- nData

  list("starts" = starts,
       "ends" = ends)
}

#' @title prep_cluster
#'
#' @description Prepares cluster for parallel computing using foreach
#' @details Prepares cluster for parallel computing using foreach
#' @param len Number of elemnts in the iterable list
#' @param n_cores Number of cores to use
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel clusterEvalQ makeCluster
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return opts and cl for foreach
#'
prep_cluster <- function(len,
                         n_cores = 2
){
  opts <- list()
  pb <- txtProgressBar(min = 0, max = len, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, .libPaths())
  doSNOW::registerDoSNOW(cl)
  list("opts" = opts, "cl" = cl)
}

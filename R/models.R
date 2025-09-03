#' @title NB_footprintings
#'
#' @description calculate footprintings by negative binomial method,
#' @details calculate pvalue and effect sizes of footprintings by negative binomial method,
#' @param Tn5Insertion a vector of Tn5 insertion for consecutive positions in the peak, matched with pred_bias, positions
#' @param pred_bias a vector of predicted Tn5 bias for consecutive positions in the peak, matched with Tn5Insertion, positions
#' @param positions a vector of consecutive positions in chr, matched with Tn5Insertion, pred_bias
#' @param center_sizes a vector of radium sizes of footprint
#' @param flank_step the step size of flanking sizes seq
#' @param p.adjust.method BH or BY, pval adjusted method
#' @param wide_flank if the flanking is 2 times of not
#' @param nCore Number of cores to use
#' @import pbmcapply
#' @return a list of adjusted pval and effect size matrix, where each row is a base pair position and each column is a center window.
#' @export
#'
NB_footprintings = function(Tn5Insertion,
                            pred_bias,
                            positions = 1:length(Tn5Insertion),
                            center_sizes = 2:100,
                            flank_window_range = 40,
                            flank_step = 4,
                            p.adjust.method = 'BH',
                            wide_flank = TRUE,
                            nCores = 2) {
  start = min(positions)
  end = max(positions)

  ### get the difference of flank windows and center size
  flank_windows = seq(-as.integer(flank_window_range/2),
                      as.integer(flank_window_range/2),
                      flank_step)

  bardata = data.frame('Tn5Insertion' = as.numeric(Tn5Insertion),
                       'pred_bias' = as.numeric(pred_bias),
                       'position' = as.numeric(positions))

  ### the range of peak
  start = min(positions)
  end = max(positions)

  ### get the results for all positions
  results = pbmcmapply(function(position_index) {
    footprintingAnalysis(center_sizes = center_sizes,
                         center_position = positions[position_index],
                         start = start,
                         end = end,
                         flank_windows = flank_windows,
                         wide_flank = wide_flank,
                         bardata = bardata)
  },
  1:length(positions),
  mc.cores = nCores)

  ### separate effect size and pval results
  p_value_matrix = t(results[1:(dim(results)[1]/2), ])
  effect_size_matrix = t(results[(1+dim(results)[1]/2):dim(results)[1], ])

  ### adjust the pval by BY on the whole regions if H0 is no footprintings in the peak
  ### in this case the pval relationship is complex and BY does not assume any relationship
  if (p.adjust.method == 'BY') {
    p_value_matrix_1 = p.adjust(p_value_matrix, method = "BY")
    p_value_matrix_1 = matrix(p_value_matrix_1,
                              nrow = dim(p_value_matrix)[1],
                              ncol = dim(p_value_matrix)[2])
  } else if (p.adjust.method == 'BH') {
    ### adjust the pval by BH on the base pair if H0 is no footprintings in one specific base pair
    ### in this case the pval are positivelty correlated and we can use BH here
    p_value_matrix_1 = t(apply(p_value_matrix, 1, p.adjust, method = "BH"))
  } else {
    stop("invalid fdr method")
  }

  rownames(p_value_matrix_1) = format(as.numeric(start:end), scientific = FALSE, trim = TRUE)
  colnames(p_value_matrix_1) = format(as.numeric(2*center_sizes), scientific = FALSE, trim = TRUE)

  ### smooth effect_size_matrix (considering the neighbourhood)
  # effect_size_matrix = apply(effect_size_matrix, 2, caTools::runmax, 5)
  # effect_size_matrix = apply(effect_size_matrix, 2, pracma::conv, 5)
  # effect_size_matrix = effect_size_matrix/(2*5)

  rownames(effect_size_matrix) = format(as.numeric(start:end), scientific = FALSE, trim = TRUE)
  colnames(effect_size_matrix) = format(as.numeric(2*center_sizes), scientific = FALSE, trim = TRUE)

  result = list()
  result[['pval']] = p_value_matrix_1
  result[['effect_size']] = effect_size_matrix
  return(result)
}

#' @title remove_intersection
#'
#' @description remove regions that overlap with the reported binding_sites
#' @details remove regions that overlap with the reported binding_sites
#' @param multiscale_data a vector of Tn5 insertion for consecutive positions in the peak, matched with pred_bias, positions
#' @param binding_sites a vector of predicted Tn5 bias for consecutive positions in the peak, matched with Tn5Insertion, positions
#' @return regions with pvalue non-overlapped with reported binding_sites
#'
remove_intersection = function(multiscale_data, binding_sites){
  if (!is.null(binding_sites)){
    return (all(multiscale_data[2] - multiscale_data[1]/2 >= binding_sites$position +
                  binding_sites$width/2|
                  multiscale_data[2] + multiscale_data[1]/2 <= binding_sites$position -
                  binding_sites$width/2)
    )
  } else {
    return(TRUE)
  }
}

#' @title binding_sites_pval
#'
#' @description searching binding sites by pvalue
#' @details searching binding sites by pvalue
#' @param p_value_matrix matrix of adjusted pvalue, where each row is a base pair position and each column is a center window.
#' @param smoothed T/F, if the pvalue needed to be smoothed within neighbourhood
#' @param pval_threshold fdr of reported binding sites
#' @param width_threshold minimum width reported binding sites
#' @import tibble
#' @importFrom tidyr pivot_longer
#' @return a list of non overlapping footprints with maximum log pvalues, rank by log pvalues
#' @export
#'
binding_sites_pval = function(p_value_matrix,
                              smoothed = FALSE,
                              pval_threshold = 0.1,
                              width_threshold = 10){
  ### log 10 pval and smooth it
  log_p_value = -log10(p_value_matrix)
  log_p_value = as.matrix(log_p_value)
  if (smoothed){
    log_p_value = apply(log_p_value, 2, caTools::runmax, 5)
    log_p_value = apply(log_p_value, 2, pracma::conv, 5)
    log_p_value = log_p_value/(2*5)
  }

  start = min(as.numeric(rownames(p_value_matrix)))
  end = max(as.numeric(rownames(p_value_matrix)))

  rownames(log_p_value) = rownames(p_value_matrix)
  log_p_value <- t(log_p_value) %>%
    as.data.frame() %>%
    rownames_to_column("width") %>%
    pivot_longer(-width, names_to = "position", values_to = "minus_log_p_value")

  log_p_value$width = as.numeric(log_p_value$width)
  log_p_value$position = as.numeric(log_p_value$position)

  ### find the candidates binding sites that pass the pval and width threshold
  multiscale_data = log_p_value %>%
    subset(minus_log_p_value >= -log10(pval_threshold) &
             width >= width_threshold &
             position >= start + width/2 + 50 &
             position <= end - width/2 - 50)
  # multiscale_data = multiscale_data[sample(1:dim(multiscale_data)[1],
  #                                          dim(multiscale_data)[1]), ]
  binding_sites = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(binding_sites) = c("width", "position", "-log_p_value")

  ### searching for max -log10 pval and remove the regions that overlap with it
  ### search by forloop
  while (dim(multiscale_data)[1] > 0) {
    candidate = multiscale_data[which.max(multiscale_data$minus_log_p_value), ]
    binding_sites = rbind(binding_sites, candidate)
    multiscale_data = multiscale_data %>%
      subset(apply(multiscale_data, 1, remove_intersection, binding_sites = binding_sites))
  }

  return(binding_sites)
}

#' @title binding_sites
#'
#' @description searching binding sites by pval
#' @details the positions and sizes of reported footprints are decided by effect size
#' @param p_value_matrix matrix of adjusted pvalues, where each row is a base pair position and each column is a center window.
#' @param effect_size_matrix matrix of effect sizes, where each row is a base pair position and each column is a center window.
#' @param pval_threshold fdr of reported binding sites
#' @param width_threshold minimum width reported binding sites
#' @import tibble
#' @importFrom caTools runmax
#' @importFrom pracma conv
#' @importFrom tidyr pivot_longer
#' @return a list of non overlapping footprints with maximum pval, but the positions and sizes of reported footprints are decided by effect size
#' @export
#'
binding_sites = function(p_value_matrix,
                         effect_size_matrix,
                         pval_threshold = 0.1,
                         width_threshold = 10){
  binding_site_candidates = binding_sites_pval(p_value_matrix,
                                               smoothed = FALSE,
                                               pval_threshold,
                                               width_threshold)
  if (dim(binding_site_candidates)[1] == 0){
    return (binding_site_candidates)
  }
  # effect_size_matrix[p_value_matrix > p_threshold &
  #                      effect_size_matrix > 0] = 0
  effect_size_matrix = apply(effect_size_matrix, 2, caTools::runmax, 5)
  effect_size_matrix = apply(effect_size_matrix, 2, pracma::conv, 5)
  effect_size_matrix = effect_size_matrix/(2*5)
  rownames(effect_size_matrix) = rownames(p_value_matrix)
  colnames(effect_size_matrix) = colnames(p_value_matrix)

  start = min(as.numeric(rownames(p_value_matrix)))
  end = max(as.numeric(rownames(p_value_matrix)))

  log_p_value = -log10(p_value_matrix)
  log_p_value = as.matrix(log_p_value)

  binding_site_candidates$effect_size = NA
  binding_site_candidates$position_effect_size = NA
  binding_site_candidates$width_effect_size = NA
  binding_site_candidates$pval_effect_size = NA

  for (j in 1:dim(binding_site_candidates)[1]){
    ranges = (binding_site_candidates$position[j] - binding_site_candidates$width[j]/2):(binding_site_candidates$position[j] + binding_site_candidates$width[j]/2)
    effect_size_matrix_range = effect_size_matrix[as.character(ranges),
                                                  as.character(2*(max(5, floor(binding_site_candidates$width[j]/3)):min(100, binding_site_candidates$width[j])))]
    binding_site_candidates$effect_size[j] = max(effect_size_matrix_range)
    binding_site_candidates[j, 5:6] = cbind(rownames(effect_size_matrix_range)[which(effect_size_matrix_range == max(effect_size_matrix_range),
                                                                                     arr.ind = T)[1, 1]],
                                            colnames(effect_size_matrix_range)[which(effect_size_matrix_range == max(effect_size_matrix_range),
                                                                                     arr.ind = T)[1, 2]])
    binding_site_candidates$pval_effect_size[j] = log_p_value[binding_site_candidates$position_effect_size[j],
                                                              binding_site_candidates$width_effect_size[j]]
  }
  return(binding_site_candidates)
}

#' @title fisher_comparison
#'
#' @description compare binding events between different conditions with fisher exact test
#' @details keep for legacy
#' @return a list of footprint binding sites with pval indicating significance of differential bindings between two conditions
fisher_comparison = function(binding_site,
                             Tn5Insertion_1,
                             Tn5Insertion_2,
                             pred_bias,
                             positions,
                             flank_window_range = 40,
                             flank_step = 2,
                             alternative = 'two.sided'){
  if (!alternative %in% c('two.sided', 'greater', 'less')){
    stop()
  }
  start = min(positions)
  end = max(positions)

  center_size = binding_site$width/2
  center_position = binding_site$position
  bardata_1 =  data.frame('Tn5Insertion' = Tn5Insertion_1,
                          'pred_bias' = pred_bias,
                          'position' = positions)
  bardata_2 =  data.frame('Tn5Insertion' = Tn5Insertion_2,
                          'pred_bias' = pred_bias,
                          'position' = positions)

  center_region = (center_position - center_size):(center_position + center_size) ### get the center region
  center_reads_1 = sum(bardata_1[bardata_1$position %in% center_region, 'Tn5Insertion'])
  center_reads_2 = sum(bardata_2[bardata_2$position %in% center_region, 'Tn5Insertion'])
  center_pred_bias = sum(bardata_1[bardata_1$position %in% center_region,  "pred_bias"]) ### get the sum of predicted bias in center region

  flank_windows = seq(-as.integer(flank_window_range/2),
                      as.integer(flank_window_range/2),
                      flank_step)

  flank_sizes = unique(pmin(pmax(center_size + flank_windows, 1), 100)) ### flank_sizes = [center_size - 40, center_size + 40]
  result = sapply(flank_sizes, function(flank_size) { ### for each flank_size
    left_region = max(start, center_position - center_size - 2*flank_size):max(start, center_position - center_size - 1)
    right_region = min(end, center_position + center_size + 1):min(end, center_position + center_size + 2*flank_size)
    left_region = setdiff(left_region, center_position)
    right_region = setdiff(right_region, center_position)
    left_index = which(bardata_1$position %in% left_region)
    right_index = which(bardata_1$position %in% right_region)
    left_reads_1 = sum(bardata_1[left_index, "Tn5Insertion"])
    right_reads_1 = sum(bardata_1[right_index, "Tn5Insertion"])
    left_reads_2 = sum(bardata_2[left_index, "Tn5Insertion"])
    right_reads_2 = sum(bardata_2[right_index, "Tn5Insertion"])
    left_pred_bias = sum(bardata_1[left_index, "pred_bias"])
    right_pred_bias = sum(bardata_1[right_index, "pred_bias"])
    ### if left or right reads = 0, make it no footprintings
    if (left_reads_1 + left_reads_2 == 0 | right_reads_1 + right_reads_2 == 0){
      pval = 1
    } else {
      ### separate test for each side
      pval_left = pnbinom(center_reads_1 + center_reads_2, size = left_reads_1 + left_reads_2,
                          prob = left_pred_bias/(left_pred_bias + center_pred_bias))
      pval_right = pnbinom(center_reads_1 + center_reads_2, size = right_reads_1 + right_reads_2,
                           prob = right_pred_bias/(right_pred_bias + center_pred_bias))
      ### conservative pval and effect size
      pval = max(pval_left, pval_right, 1e-50)
    }
    pval
  })

  if (sum(result <= 0.1) > 0){
    flank_sizes = flank_sizes[result <= 0.1]
  }

  results = sapply(flank_sizes, function(flank_size){
    flank_region_left = max(start, center_position - center_size - flank_size):max(start, center_position - center_size - 1)
    flank_region_right = min(end, center_position + center_size + 1):min(end, center_position + center_size + flank_size)

    flank_reads_left_1 = sum(bardata_1[bardata_1$position %in% flank_region_left, 'Tn5Insertion'])
    flank_reads_left_2 = sum(bardata_2[bardata_2$position %in% flank_region_left, 'Tn5Insertion'])
    flank_reads_right_1 = sum(bardata_1[bardata_1$position %in% flank_region_right, 'Tn5Insertion'])
    flank_reads_right_2 = sum(bardata_2[bardata_2$position %in% flank_region_right, 'Tn5Insertion'])

    pvals = min(fisher.test(matrix(c(center_reads_1, center_reads_2,
                                     flank_reads_left_1, flank_reads_left_2),
                                   nrow = 2), alternative = alternative)$p.value,
                fisher.test(matrix(c(center_reads_1, center_reads_2,
                                     flank_reads_right_1, flank_reads_right_2),
                                   nrow = 2), alternative = alternative)$p.value)

    if (fisher.test(matrix(c(center_reads_1, center_reads_2,
                             flank_reads_left_1, flank_reads_left_2),
                           nrow = 2), alternative = alternative)$p.value < fisher.test(matrix(c(center_reads_1, center_reads_2,
                                                                                                flank_reads_right_1, flank_reads_right_2),
                                                                                              nrow = 2), alternative = alternative)$p.value) {
      statistic = fisher.test(matrix(c(center_reads_1, center_reads_2,
                                       flank_reads_left_1, flank_reads_left_2),
                                     nrow = 2), alternative = alternative)$estimate
    } else {
      statistic = fisher.test(matrix(c(center_reads_1, center_reads_2,
                                       flank_reads_right_1, flank_reads_right_2),
                                     nrow = 2), alternative = alternative)$estimate
    }
    c(pvals, statistic)
  })
  #  weights = exp(-abs(flank_sizes - center_size))/sum(exp(- abs(flank_sizes - center_size)))
  pvals = results[1, ]
  statistics = results[2, ]
  pvals[pvals >= 1] = 1 - 1/length(pvals)
  pval = ACAT::ACAT(pvals)
  statistics = statistics[is.finite(statistics)]
  list('pval' = pval,
       'odd_ratio' = mean(statistics))
}

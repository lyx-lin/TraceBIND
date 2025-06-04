library(patchwork)
library(tidyverse)
### function for plot multiscale -log10 pval
plot_multiscale_logp = function(p_value_matrix, 
                                chr = NULL, 
                                valid_data = NULL, 
                                scale = FALSE, 
                                smooth = T){
  ### -log10 pval and smooth it
  log_p_value = -log10(p_value_matrix)
  log_p_value = as.matrix(log_p_value)
  if (smooth){
    log_p_value = apply(log_p_value, 2, caTools::runmax, 5)
    log_p_value = apply(log_p_value, 2, pracma::conv, 5)
    log_p_value = log_p_value/(2*5)
  }
  
  start = min(as.numeric(rownames(p_value_matrix)))
  end = max(as.numeric(rownames(p_value_matrix)))
  
  rownames(log_p_value) = rownames(p_value_matrix)
  colnames(log_p_value) = colnames(p_value_matrix)
  
  max_center_size = max(as.numeric(colnames(log_p_value)))
  min_center_size = min(as.numeric(colnames(log_p_value)))
  
  ### transfer the matrix into data for ggplots
  log_p_value <- t(log_p_value) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("width") %>%
    tidyr::pivot_longer(-width, names_to = "position", values_to = "p_value")
  
  log_p_value$width = as.numeric(log_p_value$width)
  log_p_value$position = as.numeric(log_p_value$position)
  
  if (scale){
    log_p_value$p_value[log_p_value$p_value > 2] = 2
    if (is.null(valid_data)){
      p_value_plot = ggplot(log_p_value, aes(position, width)) +
        geom_tile(aes(fill = p_value), alpha = 0.95) +
        xlab(paste0(chr, " center position(bp)")) +
        ylab("center size") + 
        # geom_hline(yintercept = 147, linetype='dotted', 
        #            color = 'red') + 
        # annotate("text", x = start + 50, 
        #          y=147.5, label="Nucleosome size") + 
        scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                           limits = c(start, end)) +
        scale_y_continuous(breaks = seq(0, 200, 20), 
                           limits = c(min_center_size, max_center_size)) +
        scale_fill_gradientn(name = "-log(p_val)", 
                             colours= c("white", "#F7FBFFFF", 
                                        "#DEEBF7FF", "#C6DBEFFF", 
                                        "#9ECAE1FF", "#6BAED6FF", "#4292C6FF", 
                                        "#2171B5FF", "#08519CFF", "#08306BFF"), 
                             limits = c(-0.1, 2.1), 
                             values= c(seq(0, 1, 0.1))) +
        theme_minimal() + 
        theme(panel.grid = element_line(color="black"))
    } else {
      p_value_plot = ggplot(log_p_value, aes(position, width)) +
        geom_tile(aes(fill = p_value), alpha = 0.95) +
        geom_point(data = valid_data, 
                   aes(x = x, y = y, color = source), size = 1) +
        xlab(paste0(chr, " center position(bp)")) +
        ylab("center size") + 
        # geom_hline(yintercept = 147, linetype='dotted', 
        #            color = 'red') + 
        # annotate("text", x = start + 50, 
        #          y=147.5, label="Nucleosome size") + 
        scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                           limits = c(start, end)) +
        scale_y_continuous(breaks = seq(0, 200, 20), 
                           limits = c(min_center_size, max_center_size)) +
        scale_fill_gradientn(name = "-log(p_val)", 
                             colours= c("white", "#F7FBFFFF", 
                                        "#DEEBF7FF", "#C6DBEFFF", 
                                        "#9ECAE1FF", "#6BAED6FF", "#4292C6FF", 
                                        "#2171B5FF", "#08519CFF", "#08306BFF"), 
                             limits = c(-0.1, 2.1), 
                             values = c(seq(0, max(2, log_p_value$p_value, na.rm = TRUE), length.out = 10)/max(2, log_p_value$p_value, na.rm = TRUE))) +
        theme_minimal() + 
        theme(# legend.title=element_text(size=16), 
          # legend.text=element_text(size=16), 
          panel.grid = element_line(color="black"))
    } 
  } else {
    if (is.null(valid_data)){
      p_value_plot = ggplot(log_p_value, aes(position, width)) +
        geom_tile(aes(fill = p_value), alpha = 0.95) +
        xlab(paste0(chr, " center position(bp)")) +
        ylab("center size") + 
        # geom_hline(yintercept = 147, linetype='dotted', 
        #            color = 'red') + 
        # annotate("text", x = start + 50, 
        #          y=147.5, label="Nucleosome size") + 
        scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                           limits = c(start, end)) +
        scale_y_continuous(breaks = seq(0, 200, 20), 
                           limits = c(min_center_size, max_center_size)) +
        scale_fill_gradientn(name = "-log(p_val)", 
                             colours= c("white", "#F7FBFFFF", 
                                        "#DEEBF7FF", "#C6DBEFFF", 
                                        "#9ECAE1FF", "#6BAED6FF", "#4292C6FF", 
                                        "#2171B5FF", "#08519CFF", "#08306BFF"), 
                             # limits = c(-0.1, 2.1), 
                             values= c(seq(0, 1, 0.1))) +
        theme_minimal() + 
        theme(panel.grid = element_line(color="black"))
    } else {
      p_value_plot = ggplot(log_p_value, aes(position, width)) +
        geom_tile(aes(fill = p_value), alpha = 0.95) +
        geom_point(data = valid_data, 
                   aes(x = x, y = y, color = source), size = 1) +
        xlab(paste0(chr, " center position(bp)")) +
        ylab("center size") + 
        # geom_hline(yintercept = 147, linetype='dotted', 
        #            color = 'red') + 
        # annotate("text", x = start + 50, 
        #          y=147.5, label="Nucleosome size") + 
        scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                           limits = c(start, end)) +
        scale_y_continuous(breaks = seq(0, 200, 20), 
                           limits = c(min_center_size, max_center_size)) +
        scale_fill_gradientn(name = "-log(p_val)", 
                             colours= c("white", "#F7FBFFFF", 
                                        "#DEEBF7FF", "#C6DBEFFF", 
                                        "#9ECAE1FF", "#6BAED6FF", "#4292C6FF", 
                                        "#2171B5FF", "#08519CFF", "#08306BFF"), 
                             # limits = c(-0.1, 2.1), 
                             values = c(seq(0, max(2, log_p_value$p_value, na.rm = TRUE), length.out = 10)/max(2, log_p_value$p_value, na.rm = TRUE))) +
        theme_minimal() + 
        theme(# legend.title=element_text(size=16), 
          # legend.text=element_text(size=16), 
          panel.grid = element_line(color="black"))
    }
  }
  
  p_value_plot
}

plot_multiscale_effectsize = function(effect_size_matrix, 
                                      p_value_matrix, # if p_value is not significant, then we will make the positive effect size = 0
                                      valid_data = NULL, 
                                      chr = NULL, 
                                      scale = FALSE){
  start = min(as.numeric(rownames(effect_size_matrix)))
  end = max(as.numeric(rownames(effect_size_matrix)))
  log_p_value = -log10(p_value_matrix)
  log_p_value = as.matrix(log_p_value)
  log_p_value = apply(log_p_value, 2, caTools::runmax, 5)
  log_p_value = apply(log_p_value, 2, pracma::conv, 5)
  log_p_value = log_p_value/(2*5)
  # effect_size_matrix[log_p_value <= 0.1 & 
  #                      effect_size_matrix > 0] = 0
  rownames(log_p_value) = rownames(p_value_matrix)
  
  if (scale){
    effect_size_matrix[effect_size_matrix > 0 & 
                         log_p_value <= 0.3] = 0
    # effect_size_matrix[effect_size_matrix <= 0] = 0
  }
  effect_size_matrix[, 1:4] = 0
  effect_size_matrix = apply(effect_size_matrix, 2, caTools::runmax, 5)
  effect_size_matrix = apply(effect_size_matrix, 2, pracma::conv, 5)
  effect_size_matrix = effect_size_matrix/(2*5)
  rownames(effect_size_matrix) = rownames(p_value_matrix)
  log_p_value <- t(log_p_value) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("width") %>%
    tidyr::pivot_longer(-width, names_to = "position", values_to = "p_value")
  
  log_p_value$width = as.numeric(log_p_value$width)
  log_p_value$position = as.numeric(log_p_value$position)
  
  effect_size_matrix <- t(effect_size_matrix) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("width") %>%
    tidyr::pivot_longer(-width, names_to = "position", values_to = "effect_size")
  
  effect_size_matrix$width = as.numeric(effect_size_matrix$width)
  effect_size_matrix$position = as.numeric(effect_size_matrix$position)
  # effect_size_matrix[log_p_value$p_value < -log10(0.05), 'effect_size'] = 0
  if (quantile(effect_size_matrix$effect_size, 0.99) > 0){
    effect_size_matrix[effect_size_matrix$effect_size >= quantile(effect_size_matrix$effect_size, 0.99), 'effect_size'] = quantile(effect_size_matrix$effect_size, 0.99)
  }
  if (is.null(valid_data)){
    effect_size_plot = ggplot(effect_size_matrix, aes(position, width)) +
      geom_tile(aes(fill = effect_size), alpha = 0.95) +
      xlab(paste0(chr, " center position(bp)")) +
      ylab("center size") + 
      # geom_hline(yintercept = 147, linetype='dotted', 
      #            color = 'red') + 
      # annotate("text", x = start + 50, 
      #          y=147.5, label="Nucleosome size") + 
      scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                         limits = c(start, end)) +
      scale_y_continuous(breaks = seq(0, 200, 20)) +
      scale_fill_gradient2(name = "effect size", 
                           # breaks=seq(-0.5, 0.3, 0.2), 
                           # limits = c(-0.5, 0.3), , 
                           low = scales::muted('blue'), 
                           high = scales::muted('red')) +
      theme_minimal() + 
      theme(panel.grid = element_line(color="black"))
  } else {
    effect_size_plot = ggplot(effect_size_matrix, aes(position, width)) +
      geom_tile(aes(fill = effect_size), alpha = 0.95) +
      geom_point(data = valid_data, 
                 aes(x = x, y = y, color = source), size = 1) +
      xlab(paste0(chr, " center position(bp)")) +
      ylab("center size") + 
      # geom_hline(yintercept = 147, linetype='dotted', 
      #            color = 'red') + 
      # annotate("text", x = start + 50, 
      #          y=147.5, label="Nucleosome size") + 
      scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]), 
                         limits = c(start, end)) +
      scale_y_continuous(breaks = seq(0, 200, 20)) +
      scale_fill_gradient2(name = "effect size",
                           # breaks=seq(-0.5, 0.3, 0.2), 
                           # limits = c(-0.5, 0.3), 
                           low = scales::muted('blue'), 
                           high = scales::muted('red')) +
      theme_minimal() + 
      theme(panel.grid = element_line(color="black"))
  }
  
  effect_size_plot
}

plot_insertions = function(Tn5Insertion, 
                           positions = 1:length(Tn5Insertion), 
                           chr = NULL){
  start = min(positions)
  end = max(positions)
  bardata = data.frame(position = positions, 
                       Tn5Insertion = Tn5Insertion)
  Tn5 = ggplot(bardata) + 
    geom_rect(aes(xmin = position - 2, xmax = position + 2, 
                  ymin = 0, ymax = Tn5Insertion), 
              fill = "black") +
    xlab(paste(chr ,"position (bp)")) + 
    ylab(paste0("Tn5 Insertion")) +
    scale_x_continuous(breaks = c((start:end)[(start:end) %% 100 == 0]),
                       limits = c(start, end)) +
    theme_minimal()
}

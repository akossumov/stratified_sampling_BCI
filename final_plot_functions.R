
plot_SE_boxplots <- function(n_srswor_list, n_geo_list, n_crfFe_list, SEQ_n, H){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  convert_list_to_df <- function(lst, method_name, SEQ_n) {
    bind_rows(lapply(seq_along(lst), function(i) {
      est_vals <- lst[[i]]$est_SE
      true_val <- lst[[i]]$SE
      data.frame(
        n = SEQ_n[i],
        value = est_vals,
        method = method_name,
        true_value = true_val
      )
    }))
  }
  
  
  df_all <- bind_rows(
    convert_list_to_df(n_srswor_list, "srswor", SEQ_n),
    convert_list_to_df(n_geo_list, "geo", SEQ_n),
    convert_list_to_df(n_crfFe_list, "crfFe", SEQ_n)
  )
  
  # df_all - have columns: n, value, method, true_value
  # df_true - unique true values for all (n, method)
  df_true <- df_all %>%
    select(n, method, true_value) %>%
    distinct() %>%
    mutate(n_num = as.numeric(factor(n)))  # numeric version of n for geom_segment
  
  method_colors <- c(
    "srswor" = "blue",
    "crfFe" = "green",
    "geo" = "red"
  )
  
  method_true_colors <- c(
    "srswor" = "darkblue",
    "crfFe" = "darkgreen",
    "geo" = "darkred"
  )
  
  ggplot_object <- ggplot(df_all, aes(x = factor(n), y = value, fill = method)) +
    
    #geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
    
    geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 0.8) +
    
    coord_cartesian(ylim = c(quantile(df_all$value, 0.001), quantile(df_all$value, 0.999))) +
    
    geom_segment(data = df_true,
                 aes(x = n_num - 0.35, xend = n_num + 0.35,
                     y = true_value, yend = true_value,
                     color = method),
                 size = 1.5,
                 inherit.aes = FALSE) +
    
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_true_colors) +
    
    labs(
      x = "Sample size (n)",
      y = "Standard error",
      fill = "Estimated std. errors (from simulations)",
      color = "True std. error",
      caption = paste0("NUMBER OF STRATA = ", H)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
    )
  #ggsave("SE_plot.pdf", width = 11, height = 10, units = "in")
  return(ggplot_object)
}



plot_desef_boxplots <- function(n_geo_list, n_crfFe_list, SEQ_n, H){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  convert_list_to_df <- function(lst, method_name, SEQ_n) {
    bind_rows(lapply(seq_along(lst), function(i) {
      est_vals <- lst[[i]]$est_desef
      true_val <- lst[[i]]$desef
      data.frame(
        n = SEQ_n[i],
        value = est_vals,
        method = method_name,
        true_value = true_val
      )
    }))
  }
  
  
  df_all <- bind_rows(
    convert_list_to_df(n_geo_list, "geo", SEQ_n),
    convert_list_to_df(n_crfFe_list, "crfFe", SEQ_n)
  )
  
  # df_all - have columns: n, value, method, true_value
  # df_true - unique true values for all (n, method)
  df_true <- df_all %>%
    select(n, method, true_value) %>%
    distinct() %>%
    mutate(n_num = as.numeric(factor(n)))  # numeric version of n for geom_segment
  
  method_colors <- c(
    "crfFe" = "green",
    "geo" = "red"
  )
  
  method_true_colors <- c(
    "crfFe" = "darkgreen",
    "geo" = "darkred"
  )
  
  ggplot_object <- ggplot(df_all, aes(x = factor(n), y = value, fill = method)) +
    
    #geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
    
    geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 0.8) +
    
    coord_cartesian(ylim = c(0, 1)) +
    
    geom_segment(data = df_true,
                 aes(x = n_num - 0.35, xend = n_num + 0.35,
                     y = true_value, yend = true_value,
                     color = method),
                 size = 1.5,
                 inherit.aes = FALSE) +
    
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_true_colors) +
    
    labs(
      x = "Sample size (n)",
      y = "Design effect",
      fill = "Estimated design effects (from simulations)",
      color = "True design effect",
      caption = paste0("NUMBER OF STRATA = ", H)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
    )
  #ggsave("SE_plot.pdf", width = 11, height = 10, units = "in")
  return(ggplot_object)
}


plot_mz_boxplots <- function(n_srswor_list, n_geo_list, n_crfFe_list, SEQ_n, H, mz_true){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  convert_list_to_df <- function(lst, method_name, SEQ_n) {
    bind_rows(lapply(seq_along(lst), function(i) {
      est_vals <- lst[[i]]$mz
      data.frame(
        n = SEQ_n[i],
        value = est_vals,
        method = method_name
      )
    }))
  }
  
  
  df_all <- bind_rows(
    convert_list_to_df(n_srswor_list, "srswor", SEQ_n),
    convert_list_to_df(n_geo_list, "geo", SEQ_n),
    convert_list_to_df(n_crfFe_list, "crfFe", SEQ_n)
  )
  
  
  method_colors <- c(
    "srswor" = "blue",
    "crfFe" = "green",
    "geo" = "red"
  )
  
  
  ggplot_object <- ggplot(df_all, aes(x = factor(n), y = value, fill = method)) +
    
    #geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
    
    geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 0.8) +
    
    coord_cartesian(ylim = c(quantile(df_all$value, 0.001), quantile(df_all$value, 0.999))) +
    
    geom_hline(aes(yintercept = mz_true, color = "True population mean"), linetype = "solid", linewidth = 1) +
    
    scale_fill_manual(values = method_colors) +
    scale_color_manual(
      name = "",  # name of the legend
      values = c("True population mean" = "yellow")
    ) +
    
    labs(
      x = "Sample size (n)",
      y = "Population mean",
      fill = "Estimated pop. means (from simulations)",
      caption = paste0("NUMBER OF STRATA = ", H)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
    )
  #ggsave("SE_plot.pdf", width = 11, height = 10, units = "in")
  return(ggplot_object)
}
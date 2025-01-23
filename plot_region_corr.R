plot_region_corr <- function(data, scale_name, region_name, legend_values, col_lh, col_rh, scale_name_text) {
  ggplot(data %>% filter(scale == scale_name, region == region_name)) +
    aes(x = scale_value, y = region_value, color = hemi) +
    geom_point(alpha = 0.975) +
    geom_smooth(se = FALSE, method = "lm", show.legend = FALSE)+
    theme_bw() +
    theme(
      legend.position = "top"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5))) +
    ylab("% of change") +
    xlab(scale_name_text) +
    scale_fill_manual(values = c(col_lh, col_rh), name = region_name,
                      labels = c(paste0("lh \n", legend_values[1]), 
                                 paste0("rh \n", legend_values[2]))) +
    scale_color_manual(values = c(col_lh, col_rh), name = region_name,
                       labels = c(paste0("lh \n", legend_values[1]), 
                                  paste0("rh \n", legend_values[2])))
  
}

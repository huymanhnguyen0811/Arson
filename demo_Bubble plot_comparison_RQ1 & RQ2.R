library(cowplot)

bubble_plot_rq1 <- function(df) {
  df_gas <- df %>%
    filter(., fuel_type == "Gas") %>%
    group_by(collapsed_compound) %>%
    summarise(rt1 = mean(RT1), 
              rt2 = mean(RT2),
              area = median(Percent_Area))
  
  df_diesel <- df %>%
    filter(., fuel_type == "Diesel") %>%
    group_by(collapsed_compound) %>%
    summarise(rt1 = mean(RT1), 
              rt2 = mean(RT2),
              area = median(Percent_Area))
  
  diesel_bubble <- ggplot(data = df_diesel, aes(x = rt1, y = rt2, size = area)) +
    geom_point(pch = 21) + 
    scale_size(limits = c(min(min(df_gas$area), min(df_diesel$area)), 
                          max(max(df_gas$area), max(df_diesel$area)))) +
    scale_alpha(limits = c(min(min(df_gas$area), min(df_diesel$area)), 
                           max(max(df_gas$area), max(df_diesel$area)))) +
    theme(legend.position = "none") + 
    ggtitle("Diesel")
  
  gas_bubble <- ggplot(data = df_gas, aes(x = rt1, y = rt2, size = area)) +
    geom_point(pch = 21) + 
    scale_size(limits = c(min(min(df_gas$area), min(df_diesel$area)), 
                          max(max(df_gas$area), max(df_diesel$area)))) +
    scale_alpha(limits = c(min(min(df_gas$area), min(df_diesel$area)), 
                           max(max(df_gas$area), max(df_diesel$area)))) +
    theme(legend.position = "none") + 
    ggtitle("Gas")
  
  legend <- cowplot::get_legend(diesel_bubble +
                                  theme(legend.position = "bottom"))
  
  grid.arrange(grobs = list(diesel_bubble, gas_bubble), ncol = 2, bottom = legend)
}

# Research Question 1: Gas vs. Diesel -------------------------------------------------------
rq1_statsPCA_inputdf <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% rownames(cat_5))
  
rq1_significant_wilcox_alpha0.1_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% toploadings_rq1) # alpha0.05$collapsed_compound

rq1_nonsignificant_wilcox_alpha0.1_df <- setdiff(rq1_statsPCA_inputdf, rq1_significant_wilcox_alpha0.1_df)

ASTM_rq1_cat5_wilcox_alpha0.1_names <- unique((alpha0.1 %>%
                                                        filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcoxon_stats_alpha0.1))$collapsed_compound)

ASTM_rq1_cat5_wilcox_alpha0.1_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcox_alpha0.1_names)

beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df

# Scatter plot with hollow circular points Gas vs. Diesel - compounds that significant with Wilcoxon alpha0.1
bubble_plot_rq1(rq1_significant_wilcox_alpha0.1_df)

# Scatter plot with hollow circular points Gas vs. Diesel - compounds that NOT significant with Wilcoxon alpha0.1
bubble_plot_rq1(rq1_nonsignificant_wilcox_alpha0.1_df) 

# two bubble plots supposed to look the same but didn't... sth wrong with Wilcoxon tests!!!
# ->> Investigate with histograms on values of cat_5 df -----------------------------------
# Window 1: 10<rt1<17 & 4<rt2<4.5
window1 <- rq1_nonsignificant_wilcox_alpha0.1_df %>%
  filter(., 10 < RT1 & RT1 < 17 & 4 < RT2 & RT2 < 4.5)

window1_plotlist <- list()
i <- 1
for (comp in unique(window1$collapsed_compound)) {
  window1_plotlist[[i]] <- ggplot(data = window1 %>% filter(., collapsed_compound == comp),
         aes(x = Area, fill = fuel_type)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity") + 
    ggtitle(comp) +
    theme(legend.position = "none")
  i <- i + 1
}

legend <- cowplot::get_legend(window1_plotlist[[1]] + theme(legend.position = "bottom")) 
grid.arrange(grobs = window1_plotlist, ncol = 4, bottom = legend)
  
# Window 2: rt1<10 & 3<rt2<3.7
window2 <- rq1_nonsignificant_wilcox_alpha0.1_df %>%
  filter(., RT1 < 10 & 3 < RT2 & RT2 < 3.7)

window2_plotlist <- list()
i <- 1
for (comp in unique(window2$collapsed_compound)) {
  window2_plotlist[[i]] <- ggplot(data = window2 %>% filter(., collapsed_compound == comp),
                                  aes(x = Area, fill = fuel_type)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity") + 
    ggtitle(comp) +
    theme(legend.position = "none")
  i <- i + 1
}

legend <- cowplot::get_legend(window2_plotlist[[1]] + theme(legend.position = "bottom")) 
grid.arrange(grobs = window2_plotlist, ncol = 6, bottom = legend)

# Scatter plot with hollow circular points Gas vs. Diesel - ASTM compounds that significant with Wilcoxon alpha0.1

bubble_plot_rq1(ASTM_rq1_cat5_wilcox_alpha0.1_df)

# Scatter plot with hollow circular points Gas vs. Diesel - BEYOND ASTM compounds that significant with Wilcoxon alpha0.1

bubble_plot_rq1(beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df)

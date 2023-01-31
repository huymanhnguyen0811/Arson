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
  
rq1_significant_wilcox_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% alpha0.05$collapsed_compound)

rq1_toploading_pca <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% toploadings_rq1)

rq1_nonsignificant_wilcox_df <- setdiff(rq1_statsPCA_inputdf, rq1_significant_wilcox_df)

# ASTM_rq1_cat5_wilcox_alpha0.1_names <- unique((alpha0.1 %>%
#                                                         filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcoxon_stats_alpha0.1))$collapsed_compound)
# 
# ASTM_rq1_cat5_wilcox_alpha0.1_df <- shared_comp_normalized %>%
#   filter(., fuel_type %in% c("Gas", "Diesel")) %>%
#   filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcox_alpha0.1_names)
# 
# beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df

# Gas vs. Diesel - compounds that significant with Wilcoxon alpha0.1
bubble_plot_rq1(rq1_significant_wilcox_df)

# Gas vs. Diesel - PCA top 100 loadings
bubble_plot_rq1(rq1_toploading_pca)

# Gas vs. Diesel - compounds that NOT significant with Wilcoxon alpha0.1
bubble_plot_rq1(rq1_nonsignificant_wilcox_df) 

# two bubble plots supposed to look the same but didn't... sth wrong with Wilcoxon tests!!!
# ->> Investigate with histograms on values of cat_5 df
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

# Gas vs. Diesel - ASTM compounds that significant with Wilcoxon alpha0.1
bubble_plot_rq1(ASTM_rq1_cat5_wilcox_alpha0.1_df)

# Gas vs. Diesel - BEYOND ASTM compounds that significant with Wilcoxon alpha0.1
bubble_plot_rq1(beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df)

# Research Question 2: Gas stations 1,3,7 vs. 5,7,9 -------------------------------------------------------
rq2_statsPCA_inputdf <- shared_comp_normalized %>% 
  filter(., fuel_type %in% "Gas") %>%
  filter(., collapsed_compound %in% colnames(rq2_cat2_stats))

rq2_significant_wilcox_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% "Gas") %>%
  filter(., collapsed_compound %in% rq2_cat2_alpha0.05$collapsed_compound)

# rq2_toploading_pca

rq2_nonsignificant_wilcox_df <- setdiff(rq2_statsPCA_inputdf, rq2_significant_wilcox_df)

df_group1 <- rq2_significant_wilcox_df %>%
  filter(., gas_station %in% c("Station_5", "Station_7", "Station_9")) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2),
            area = median(Percent_Area))

df_group2 <- rq2_significant_wilcox_df %>%
  filter(., gas_station %in% c("Station_1", "Station_3", "Station_8")) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2),
            area = median(Percent_Area))

group1_bubble <- ggplot(data = df_group1, aes(x = rt1, y = rt2, size = area)) +
  geom_point(pch = 21) + 
  scale_size(limits = c(min(min(df_group1$area), min(df_group2$area)), 
                        max(max(df_group1$area), max(df_group2$area)))) +
  scale_alpha(limits = c(min(min(df_group1$area), min(df_group2$area)), 
                         max(max(df_group1$area), max(df_group2$area)))) +
  theme(legend.position = "none") + 
  ggtitle("Station 5,7,9")

group2_bubble <- ggplot(data = df_group2, aes(x = rt1, y = rt2, size = area)) +
  geom_point(pch = 21) + 
  scale_size(limits = c(min(min(df_group1$area), min(df_group2$area)), 
                        max(max(df_group1$area), max(df_group2$area)))) +
  scale_alpha(limits = c(min(min(df_group1$area), min(df_group2$area)), 
                         max(max(df_group1$area), max(df_group2$area)))) +
  theme(legend.position = "none") + 
  ggtitle("Station 1,3,8")

legend <- cowplot::get_legend(group1_bubble +
                                theme(legend.position = "bottom"))

grid.arrange(grobs = list(group1_bubble, group2_bubble), ncol = 2, bottom = legend)

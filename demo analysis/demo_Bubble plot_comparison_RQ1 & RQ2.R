library(cowplot)
library(ggpubr)

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
# rq1_statsPCA_inputdf <- shared_comp_normalized_rt10.2 %>%
#   filter(., fuel_type %in% c("Gas", "Diesel")) %>%
#   filter(., collapsed_compound %in% rownames(df_pca_rq1_rt10.2[[1]]))
  
# Gas significant ------------------
rq1_gas_significant <- cat_5[,1:21] %>% 
  mutate(area = rowMeans(select(., 1:21))) %>% 
  rownames_to_column(., var = "collapsed_compound") %>%
  relocate(area, .after = collapsed_compound) %>%
  select(., 1:2) %>%
  filter(., collapsed_compound %in% rq1_p0.1$collapsed_compound)

rq1_gas_significant_rt <- shared_comp_rt10.1 %>%
  filter(., fuel_type == "Gas") %>%
  filter(., collapsed_compound %in% rq1_p0.1$collapsed_compound) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2))

rq1_gas_significant <- dplyr::right_join(rq1_gas_significant_rt, rq1_gas_significant, by = "collapsed_compound")

# Diesel significant -----------
rq1_diesel_significant <- cat_5[,22:25] %>% 
  mutate(area = rowMeans(select(., 1:4))) %>% 
  rownames_to_column(., var = "collapsed_compound") %>%
  relocate(area, .after = collapsed_compound) %>%
  select(., 1:2) %>%
  filter(., collapsed_compound %in% rq1_p0.1$collapsed_compound)

rq1_diesel_significant_rt <- shared_comp_rt10.1 %>%
  filter(., fuel_type == "Diesel") %>%
  filter(., collapsed_compound %in% rq1_p0.1$collapsed_compound) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2))

rq1_diesel_significant <- dplyr::right_join(rq1_diesel_significant_rt, rq1_diesel_significant, by = "collapsed_compound")

# Top PCA loadings
# rq1_toploading_pca <- shared_comp_normalized %>%
#   filter(., fuel_type %in% c("Gas", "Diesel")) %>%
#   filter(., collapsed_compound %in% toploadings_rq1)

###### Gas Non-significant ----------------
# rq1_gas_nonsignificant <- cat_5[,1:21] %>% 
#   mutate(area = rowMeans(select(., 1:21))) %>% 
#   rownames_to_column(., var = "collapsed_compound") %>%
#   relocate(area, .after = collapsed_compound) %>%
#   select(., 1:2) %>%
#   filter(., collapsed_compound %notin% alpha0.1$collapsed_compound)
# 
# rq1_gas_nonsignificant_rt <- shared_comp_normalized_rt10.2 %>%
#   filter(., fuel_type == "Gas") %>%
#   filter(., collapsed_compound %in% rq1_gas_nonsignificant$collapsed_compound) %>%
#   group_by(collapsed_compound) %>%
#   summarise(rt1 = mean(RT1), 
#             rt2 = mean(RT2))
# 
# rq1_gas_nonsignificant <- dplyr::right_join(rq1_gas_nonsignificant_rt, rq1_gas_nonsignificant, by = "collapsed_compound")

###### Diesel Non-significant -------------
# rq1_diesel_nonsignificant <- cat_5[,22:25] %>% 
#   mutate(area = rowMeans(select(., 1:4))) %>% 
#   rownames_to_column(., var = "collapsed_compound") %>%
#   relocate(area, .after = collapsed_compound) %>%
#   select(., 1:2) %>%
#   filter(., collapsed_compound %notin% alpha0.1$collapsed_compound)
# 
# rq1_diesel_nonsignificant_rt <- shared_comp_normalized_rt10.2 %>%
#   filter(., fuel_type == "Diesel") %>%
#   filter(., collapsed_compound %in% rq1_diesel_nonsignificant$collapsed_compound) %>%
#   group_by(collapsed_compound) %>%
#   summarise(rt1 = mean(RT1), 
#             rt2 = mean(RT2))
# 
# rq1_diesel_nonsignificant <- dplyr::right_join(rq1_diesel_nonsignificant_rt, rq1_diesel_nonsignificant, by = "collapsed_compound")

# ASTM_rq1_cat5_wilcox_alpha0.1_names <- unique((alpha0.1 %>%
#                                                         filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcoxon_stats_alpha0.1))$collapsed_compound)
# 
# ASTM_rq1_cat5_wilcox_alpha0.1_df <- shared_comp_normalized %>%
#   filter(., fuel_type %in% c("Gas", "Diesel")) %>%
#   filter(., collapsed_compound %in% ASTM_rq1_cat5_wilcox_alpha0.1_names)
# 
# beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df

# Gas vs. Diesel - compounds that significant with Wilcoxon alpha0.1---------------
rq1_diesel_significant_new <- rq1_diesel_significant %>% mutate(., fuel_type = "Diesel")
rq1_gas_significant_new <- rq1_gas_significant %>% mutate(., fuel_type = "Gas")

# rq1_sig <- rbind(rq1_diesel_significant, rq1_gas_significant)

# # Add ASTM group
ASTM_group <- ASTM_list %>% filter(., !is.na(Group))

for (group in unique(ASTM_group$Group)) {
  minrt1 <- min(ASTM_group[ASTM_group$Group == group,]$RT1)
  maxrt1 <- max(ASTM_group[ASTM_group$Group == group,]$RT1)
  minrt2 <- min(ASTM_group[ASTM_group$Group == group,]$RT2)
  maxrt2 <- max(ASTM_group[ASTM_group$Group == group,]$RT2)
  idx <- which(rq1_gas_significant_new$rt1 <= maxrt1 &
                 rq1_gas_significant_new$rt1 >= minrt1 &
                 rq1_gas_significant_new$rt2 <= maxrt2 &
                 rq1_gas_significant_new$rt2 >= minrt2)
  rq1_gas_significant_new[idx, "group"] <- group
}

# Add column name
colnames(rq1_diesel_significant) <- c("collapsed_compound", "Retention time 1", "Retention time 2", 
                       "TSN-normalized area", "fuel_type") #  "Chemical Group"
colnames(rq1_gas_significant) <- c("collapsed_compound", "Retention time 1", "Retention time 2", 
                                          "TSN-normalized area", "fuel_type") #  "Chemical Group"

# MAking BUbble plot of significant RQ1 compounds for Gasoline and Diesel ===================================
  ggplot() +
    geom_point(data = rq1_diesel_significant, aes(x = `Retention time 1`, y = `Retention time 2`,
                                   size = `TSN-normalized area`), #, color = `Chemical Group`
               pch = 21) + 
  
    # stat_ellipse() +
    # facet_grid(~fuel_type) +
    scale_x_continuous(breaks = seq(from = 0, to = 40, by = 5)) +
    scale_size(limits  = c(min(min(rq1_gas_significant$`TSN-normalized area`), 
                               min(rq1_diesel_significant$`TSN-normalized area`)), 
                          max(max(rq1_gas_significant$`TSN-normalized area`), 
                              max(rq1_diesel_significant$`TSN-normalized area`))),
               range = c(2,12)) +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme_classic(base_size = 15) +
    labs(title = "Diesel") + # "Gas"
    theme(plot.title = element_text(hjust = 0.5))
# =======================================================
# plot +geom_point(data = ASTM_F001A_aligned, pch = 4,
#                  aes(x=RT1_aligned, y= RT2_aligned, colour = "red"))

# Gas vs. Diesel - PCA top 100 loadings---------
bubble_plot_rq1(rq1_toploading_pca)

# Gas vs. Diesel - compounds that NOT significant with Wilcoxon alpha0.1------------
diesel_bubble_nonsf <- ggplot(data = rq1_diesel_nonsignificant  , aes(x = rt1, y = rt2, size = area)) +
  geom_point(pch = 21) + 
  scale_size(limits = c(min(min(rq1_gas_nonsignificant $area), min(rq1_diesel_nonsignificant$area)), 
                        max(max(rq1_gas_nonsignificant $area), max(rq1_diesel_nonsignificant$area)))) +
  scale_alpha(limits = c(min(min(rq1_gas_nonsignificant $area), min(rq1_diesel_nonsignificant$area)), 
                         max(max(rq1_gas_nonsignificant $area), max(rq1_diesel_nonsignificant$area)))) +
  theme(legend.position = "none") + 
  ggtitle("Diesel_Non-significant")

gas_bubble_nonsf <- ggplot(data = rq1_gas_nonsignificant , aes(x = rt1, y = rt2, size = area)) +
  geom_point(pch = 21) + 
  scale_size(limits = c(min(min(rq1_gas_nonsignificant $area), min(rq1_diesel_nonsignificant$area)), 
                        max(max(rq1_gas_nonsignificant $area), max(rq1_diesel_nonsignificant$area)))) +
  scale_alpha(limits = c(min(min(rq1_gas_nonsignificant $area), min(rq1_diesel_nonsignificant$area)), 
                         max(max(rq1_gas_nonsignificant $area), max(rq1_diesel_nonsignificant$area)))) +
  theme(legend.position = "none") + 
  ggtitle("Gas_Non-significant")

legend <- cowplot::get_legend(diesel_bubble_nonsf +
                                theme(legend.position = "bottom"))

grid.arrange(grobs = list(diesel_bubble_nonsf, gas_bubble_nonsf), ncol = 2, bottom = legend)

# two bubble plots supposed to look the same but didn't... sth wrong with Wilcoxon tests!!!
# ->> Investigate with histograms on values of cat_5 df -------------
# Joined non-significant Gas and Diesel dfs

join_nonsf <- rbind(rq1_gas_nonsignificant %>% 
                      mutate(., fuel_type = "Gas"), 
                    rq1_diesel_nonsignificant %>% 
                      mutate(., fuel_type = "Diesel"))

# Window 1: 10<rt1<17 & 4<rt2<4.5
window1 <- join_nonsf %>%
  filter(., 10 < rt1 & rt1 < 17 & 4 < rt2 & rt2 < 4.5)

transpose_window1 <- data.table::transpose(cat_5 %>%
                                             rownames_to_column(., var = "comps") %>% 
                                             filter(., comps %in% unique(window1$collapsed_compound)) %>%
                                             column_to_rownames(., var = "comps"))

colnames(transpose_window1) <- rownames(cat_5 %>%
                                          rownames_to_column(., var = "comps") %>% 
                                          filter(., comps %in% unique(window1$collapsed_compound)) %>%
                                          column_to_rownames(., var = "comps"))
rownames(transpose_window1) <- colnames(cat_5 %>%
                                          rownames_to_column(., var = "comps") %>% 
                                          filter(., comps %in% unique(window1$collapsed_compound)) %>%
                                          column_to_rownames(., var = "comps"))

window1_plotdat <- transpose_window1 %>% 
  rownames_to_column(., var = "sample_name") %>% 
  mutate(., fuel_type = ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")) %>%
  relocate(fuel_type, .after = sample_name)

window1_plotlist <- list()
i <- 1
for (comp in colnames(window1_plotdat[, 3:ncol(window1_plotdat)])) {
  dat <- window1_plotdat %>% select(., c(comp, fuel_type))
  window1_plotlist[[i]] <- ggplot(data = dat,
         aes(x = dat[,1], fill = fuel_type)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity") + 
    ggtitle(comp) +
    theme(legend.position = "none",
          axis.title.x=element_text()) +
    xlab("Percent Area")
  i <- i + 1
}

legend <- cowplot::get_legend(window1_plotlist[[1]] + theme(legend.position = "bottom")) 
grid.arrange(grobs = window1_plotlist, ncol = 4, bottom = legend)
  
# Window 2: rt1<10 & 3<rt2<3.7
window2 <- join_nonsf %>%
  filter(., rt1 < 10 & 3 < rt2 & rt2 < 3.7)

transpose_window2 <- data.table::transpose(cat_5 %>%
                                             rownames_to_column(., var = "comps") %>% 
                                             filter(., comps %in% unique(window2$collapsed_compound)) %>%
                                             column_to_rownames(., var = "comps"))

colnames(transpose_window2) <- rownames(cat_5 %>%
                                          rownames_to_column(., var = "comps") %>% 
                                          filter(., comps %in% unique(window2$collapsed_compound)) %>%
                                          column_to_rownames(., var = "comps"))
rownames(transpose_window2) <- colnames(cat_5 %>%
                                          rownames_to_column(., var = "comps") %>% 
                                          filter(., comps %in% unique(window2$collapsed_compound)) %>%
                                          column_to_rownames(., var = "comps"))

window2_plotdat <- transpose_window2 %>% 
  rownames_to_column(., var = "sample_name") %>% 
  mutate(., fuel_type = ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")) %>%
  relocate(fuel_type, .after = sample_name)

window2_plotlist <- list()
i <- 1
for (comp in colnames(window2_plotdat[, 3:ncol(window2_plotdat)])) {
  dat <- window2_plotdat %>% select(., c(comp, fuel_type))
  window2_plotlist[[i]] <- ggplot(data = dat,
                                  aes(x = dat[,1], fill = fuel_type)) +
    geom_histogram(bins = 100, alpha = 0.5, position = "identity") + 
    ggtitle(comp) +
    theme(legend.position = "none",
          axis.title.x=element_text()) +
    xlab("Percent Area")
  i <- i + 1
}

legend <- cowplot::get_legend(window2_plotlist[[1]] + theme(legend.position = "bottom")) 
grid.arrange(grobs = window2_plotlist, ncol = 6, bottom = legend)

# Gas vs. Diesel - ASTM compounds that significant with Wilcoxon alpha0.1 ------
# bubble_plot_rq1(ASTM_rq1_cat5_wilcox_alpha0.1_df)

# Gas vs. Diesel - BEYOND ASTM compounds that significant with Wilcoxon alpha0.1 ------
# bubble_plot_rq1(beyond_ASTM_rq1_cat5_wilcox_alpha0.1_df)

# Research Question 2: Gas stations 1,3,8 vs. 5,7,9 -------------------------------------------------------
# Gas station 1,3,8 significant ------------------
df_138 <- df_stats_rq2_rt10.1[which(df_stats_rq2_rt10.1$gas_station %in% c("Station_1", "Station_2", "Station_5")), 
                                                    c(1, 3:ncol(df_stats_rq2_rt10.1))]
rownames(df_138) <- NULL
df_138 <- df_138 %>%
  column_to_rownames(., var = "sample_name")

transpose_df138 <- data.table::transpose(df_138)
rownames(transpose_df138) <- colnames(df_138)
colnames(transpose_df138) <- rownames(df_138)

rq2_gas_138 <- transpose_df138 %>% 
  mutate(area = base::rowMeans(select(., 1:12))) %>% 
  rownames_to_column(., var = "collapsed_compound") %>%
  relocate(area, .after = collapsed_compound) %>%
  select(., 1:2) %>%
  filter(., collapsed_compound %in% rq2_cat2_alpha0.1$collapsed_compound) # rq2_cat2_alpha0.1$collapsed_compound

rq2_gas_138_rt <- shared_comp_rt10.1 %>%
  filter(., fuel_type == "Gas") %>%
  filter(., collapsed_compound %in% rq2_cat2_alpha0.1$collapsed_compound) %>% # rq2_cat2_alpha0.1$collapsed_compound
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2))

rq2_gas_138_significant <- dplyr::right_join(rq2_gas_138, rq2_gas_138_rt, by = "collapsed_compound")

# Gas station 5,7,9 significant -----------
df_579 <- df_stats_rq2_rt10.1[which(df_stats_rq2_rt10.1$gas_station %in% c("Station_3", "Station_4", "Station_6")), 
                              c(1, 3:ncol(df_stats_rq2_rt10.1))]
rownames(df_579) <- NULL
df_579 <- df_579 %>%
  column_to_rownames(., var = "sample_name")

transpose_df579 <- data.table::transpose(df_579)
rownames(transpose_df579) <- colnames(df_579)
colnames(transpose_df579) <- rownames(df_579)

rq2_gas_579 <- transpose_df579 %>% 
  mutate(area = base::rowMeans(select(., 1:9))) %>% 
  rownames_to_column(., var = "collapsed_compound") %>%
  relocate(area, .after = collapsed_compound) %>%
  select(., 1:2) %>%
  filter(., collapsed_compound %in% rq2_cat2_alpha0.1$collapsed_compound) # rq2_cat2_alpha0.1$collapsed_compound

rq2_gas_579_rt <- shared_comp_rt10.1 %>%
  filter(., fuel_type == "Gas") %>%
  filter(., collapsed_compound %in% rq2_cat2_alpha0.1$collapsed_compound) %>% # rq2_cat2_alpha0.1$collapsed_compound
  group_by(collapsed_compound) %>%
  summarise(rt1 = mean(RT1), 
            rt2 = mean(RT2))

rq2_gas_579_significant <- dplyr::right_join(rq2_gas_579, rq2_gas_579_rt, by = "collapsed_compound")

# Gas 138 vs. Gas 579 - compounds that significant with Wilcoxon alpha0.1---------------
# Add column 'gas_station' to each dfs
rq2_138_significant_new <- rq2_gas_138_significant %>% mutate(., gas_station = "Gas station 1,2,5")
rq2_579_significant_new <- rq2_gas_579_significant %>% mutate(., gas_station = "Gas station 3,4,6")

# rq2_sig <- rbind(rq2_138_significant, rq2_579_significant)

# Add ASTM group
# ASTM_group <- ASTM_list %>% filter(., !is.na(Group))
# 
# for (group in unique(ASTM_group$Group)) {
#   minrt1 <- min(ASTM_group[ASTM_group$Group == group,]$RT1)
#   maxrt1 <- max(ASTM_group[ASTM_group$Group == group,]$RT1)
#   minrt2 <- min(ASTM_group[ASTM_group$Group == group,]$RT2)
#   maxrt2 <- max(ASTM_group[ASTM_group$Group == group,]$RT2)
#   idx <- which(rq2_579_significant_new$rt1 < maxrt1 &
#                  rq2_579_significant_new$rt1 > minrt1 &
#                  rq2_579_significant_new$rt2 < maxrt2 &
#                  rq2_579_significant_new$rt2 > minrt2)
#   rq2_579_significant_new[idx, "group"] <- group
# }

colnames(rq2_138_significant_new) <- c("collapsed_compound", "TSN-normalized area", "Retention time 1",
                       "Retention time 2", "gas_station") # , "Chemical Group"
colnames(rq2_579_significant_new) <- c("collapsed_compound", "TSN-normalized area", "Retention time 1",
                       "Retention time 2", "gas_station") # , "Chemical Group"

# BUUBLE PLOT
ggplot() +
  geom_point(data = rq2_138_significant_new, aes(x = `Retention time 1`, y = `Retention time 2`,
                                 size = `TSN-normalized area`), # , color = `Chemical Group`
             pch = 21) +
  # stat_ellipse() +
  # facet_grid(~gas_station) +
  scale_x_continuous(breaks = seq(from = 0, to = 40, by = 5)) +
  scale_size(limits = c(min(min(rq2_138_significant_new$`TSN-normalized area`), 
                            min(rq2_579_significant_new$`TSN-normalized area`)), 
                        max(max(rq2_138_significant_new$`TSN-normalized area`), 
                            max(rq2_579_significant_new$`TSN-normalized area`))),
             range = c(2,12)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_classic2(base_size = 15) +
  labs(title = "Gas station 1, 2, 5") + # "Gas station 3, 4, 6"
  theme(plot.title = element_text(hjust = 0.5))


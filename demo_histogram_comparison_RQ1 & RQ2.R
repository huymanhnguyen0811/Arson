# Research Question 1: Gas vs. Diesel -------------------------------------------------------
rq1_statsPCA_inputdf <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>% 
  filter(., collapsed_compound %in% rownames(cat_5))

rq1_significant_wilcox_df <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  filter(., collapsed_compound %in% alpha0.05$collapsed_compound) # 

rq1_nonsignificant_wilcox_df <- setdiff(rq1_statsPCA_inputdf, rq1_significant_wilcox_df)

######################################################

ggplot(data = rq1_significant_wilcox_df, aes(x = Percent_Area, fill = fuel_type)) +
  geom_histogram(alpha = 0.5,bins=200, position = "identity") + 
  facet_wrap(.~collapsed_compound) 

ggplot(data = rq1_significant_wilcox_df, aes(x = fuel_type, y = Percent_Area, color = fuel_type)) +
  # geom_violin(aes(fill = fuel_type), alpha = 0.05) + 
  geom_boxplot(aes(fill = fuel_type), alpha = 0.05) +
  facet_wrap(.~collapsed_compound)

ggplot(data = rq1_nonsignificant_wilcox_df, aes(x = Percent_Area, fill = fuel_type)) +
  geom_histogram(alpha = 0.5,bins=200, position = "identity") + 
  facet_wrap(.~collapsed_compound) 

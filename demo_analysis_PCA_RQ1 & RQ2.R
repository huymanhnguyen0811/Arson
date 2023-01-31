# Research Question 1 PCA -------------------------
# category 5: >=2 Gas and /or >=2 Diesel data pts
df_X_rq1 <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  select(sample_name, collapsed_compound, Percent_Area) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(sample_name, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
  relocate(`0220F001D.xlsx`, `0220F009D.xlsx`, `0220F009-2D.xlsx`, `0220F005D.xlsx`, .after = last_col()) %>%
  column_to_rownames(., var = "collapsed_compound")

# table for information ((rows are sample IDs, columns are sample information) -----------------------
metadata_X_rq1 <- data.frame(unique((shared_comp_normalized %>%
                                   filter(., fuel_type %in% c("Gas", "Diesel")))$sample_name)) 
colnames(metadata_X_rq1) <- c('sample_name')
metadata_X_rq1 <- metadata_X_rq1 %>%
  mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                            ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                   ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))) %>%
  # # Grouping samples into respective Gas stations
  mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                              ifelse(str_detect(sample_name, "F001"), "Station_1",
                                     ifelse(str_detect(sample_name, "F007"), "Station_7",
                                            ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                   ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                          ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
  column_to_rownames(., var = "sample_name")

move_to_last <- function(df, n) df[c(setdiff(seq_len(nrow(df)), n), n), ]
metadata_X_rq1 <- move_to_last(metadata_X_rq1, 5)
metadata_X_rq1 <- move_to_last(metadata_X_rq1, 2)
metadata_X_rq1 <- move_to_last(metadata_X_rq1, 7)
metadata_X_rq1 <- move_to_last(metadata_X_rq1, 21)

cat_5 <- stats_rq1[(rowSums(!is.na(stats_rq1[, 1:21])) >= 2) & 
                     (rowSums(!is.na(stats_rq1[, 22:25])) >= 2),] # -> 374 out of 4769 compounds

for (r in 1:nrow(cat_5)) { 
  cat_5[r, which(base::is.na(cat_5[r,]))] <- runif(length(which(base::is.na(cat_5[r,]))),
                                                   min = sort(shared_comp_normalized$Percent_Area)[1],
                                                   max = sort(shared_comp_normalized$Percent_Area)[2])
}

all(colnames(cat_5) == rownames(metadata_X_rq1))

# Conduct principal component analysis (PCA): ---------------------------------
p_rq1 <- pca(cat_5, metadata = metadata_X_rq1) # Classify on sample

# A bi-plot -------------
biplot(p_rq1,
       lab = row.names(p_rq1$metadata),
       colby = 'fuel_type',
       hline = 0, vline = 0,
       legendPosition = 'right', labSize = 5,
       sizeLoadingsNames = 5,
       showLoadings = TRUE,
       ntopLoadings = 10,
       pointSize = 4, 
       legendLabSize = 15,
       legendTitleSize = 16,
       legendIconSize = 6)

# Retrieve compound name of top 100 loading
loadingS_rq1_sorted <- p_rq1$loadings %>% arrange(PC1, PC2)
toploadings_rq1 <- rownames(loadingS_rq1_sorted[c(1:50, nrow(loadingS_rq1_sorted):(nrow(loadingS_rq1_sorted) - 50)),1:2])

# Cross check similarity between top 100 loading and RQ1 Wilcoxon test result ASTM_alpha0.1 - 12 compounds
idx <- which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, PASS WILCOXON TEST (data imputed with LOD), alpha threshold < 0.1 (rt1thres = 0.2)`, 
                           "x"))

name <- unlist(flatten(ASTM_list[idx, 16]))

ASTM_rq1_cat5_wilcoxon_stats_alpha0.1 <- gsub("x_", "", name)

intersect_PCA_Wilcoxon_alpha0.1 <- intersect(toploadings_rq1, ASTM_rq1_cat5_wilcoxon_stats_alpha0.1) # 15 ASTM compounds similar between top 200 loading PCA and ASTM_alpha0.1


# Pairs plot -------------
pairsplot(p_rq1,
          components = getComponents(p, c(1:5)),
          triangle = FALSE,
          trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 1.5,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'fuel_type',
          title = 'Pairs plot',
          axisLabSize = 14, plotaxes = TRUE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


# explore further the collapsed_compounds that are driving these differences along each PC. -------------
plotloadings(p_rq1,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 4)


eigencorplot(p, metavars = c('fuel_type', 'gas_station'))



# Research Question 2 PCA ################################################################################################
# Category 2: Compounds that have >=2 data points for each gas stations
# transpose the rows and columns
rq2_pca <- rq2_cat2_stats_imputed[,-2] %>% column_to_rownames(., var = "sample_name") # rq2_cat2_stats is from RQ2 of demo_analysis_stats R script
df_X_rq2 <- data.table::transpose(rq2_pca)
rownames(df_X_rq2) <- colnames(rq2_pca)
colnames(df_X_rq2) <- rownames(rq2_pca)

# table for information ((rows are sample IDs, columns are sample information) -----------------------
metadata_X_rq2 <- data.frame(unique((shared_comp_normalized %>%
                                   filter(., fuel_type %in% "Gas"))$sample_name)) 
colnames(metadata_X_rq2) <- c('sample_name')

metadata_X_rq2 <- metadata_X_rq2 %>%
  # # Grouping samples into respective Gas stations
  mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                              ifelse(str_detect(sample_name, "F001"), "Station_1",
                                     ifelse(str_detect(sample_name, "F007"), "Station_7",
                                            ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                   ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                          ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
  column_to_rownames(., var = "sample_name") 

# check that sample names match exactly between pdata and expression data 
all(colnames(df_X_rq2) == rownames(metadata_X_rq2))

# Conduct principal component analysis (PCA): ---------------------------------
p_rq2 <- pca(df_X_rq2, metadata = metadata_X_rq2)

# A bi-plot -------------
biplot(p_rq2,
       lab = row.names(p_rq2$metadata),
       colby = 'gas_station',
       hline = 0, vline = 0,
       legendPosition = 'right',labSize = 5,
       sizeLoadingsNames = 5,
       showLoadings = TRUE,
       pointSize = 4, 
       legendLabSize = 15,
       legendTitleSize = 16,
       legendIconSize = 6)

# Pairs plot -------------
pairsplot(p_rq2,
          components = getComponents(p_rq2, c(1:5)),
          triangle = FALSE,
          trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 1.5,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'gas_station',
          title = 'Pairs plot',
          axisLabSize = 14, plotaxes = TRUE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


# explore further the collapsed_compounds that are driving these differences along each PC. ------------
plotloadings(p_rq2,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 4)

eigencorplot(p_rq2, metavars = c('fuel_type', 'gas_station'))

# Retrieve compound name of top x loading
loadingS_rq2_sorted <- p_rq2$loadings %>% arrange(PC1, PC2)
toploadings_rq2_pca <- rownames(loadingS_rq2_sorted[c(1:100, nrow(loadingS_rq2_sorted):(nrow(loadingS_rq2_sorted) - 100)),1:2])

# Cross check similarity between top 100 loading and RQ1 Wilcoxon test result ASTM_alpha0.1 -  compounds
idx <- which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION, AFTER Wilcoxon test (data imputed with LOD), alpha threshold < 0.1 (rt1_thres = 0.2)`, 
                        "x"))
name <- unlist(flatten(ASTM_list[idx, 20]))

ASTM_rq2_cat2_wilcoxon_stats_alpha0.1 <- gsub("x_", "", name)

intersect_PCA_Wilcoxon_alpha0.1_rq2 <- intersect(toploadings_rq2_pca, ASTM_rq2_cat2_wilcoxon_stats_alpha0.1)


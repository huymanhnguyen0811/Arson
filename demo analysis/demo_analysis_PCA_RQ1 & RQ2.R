library(PCAtools)

# Research Question 1 PCA -----------------
# category 5 >=2 Gas and /or >=2 Diesel data pts 
df_pca_rq1 <- function(data) {
  df_X_rq1 <- data %>%
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
  
  cat_5 <- df_X_rq1[(rowSums(!is.na(df_X_rq1[, 1:21])) >= 2) & 
                      (rowSums(!is.na(df_X_rq1[, 22:25])) >= 2),] # -> 374 out of 4769 compounds
  
  for (r in 1:nrow(cat_5)) { 
    cat_5[r, which(base::is.na(cat_5[r,]))] <- runif(length(which(base::is.na(cat_5[r,]))),
                                                     min = sort(data$Percent_Area)[1],
                                                     max = sort(data$Percent_Area)[2])
  }
  
  # table for information ((rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq1 <- data.frame(colnames(cat_5)) 
  colnames(metadata_X_rq1) <- c('sample_name')
  metadata_X_rq1 <- metadata_X_rq1 %>%
    mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                              ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                     ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))) %>%
    column_to_rownames(., var = "sample_name")
  
  print(all(colnames(cat_5) == rownames(metadata_X_rq1)))
  
  return(list(cat_5 ,metadata_X_rq1))
}

df_pca_rq1_rt10.1 <- df_pca_rq1(shared_comp_rt10.1)
# df_pca_rq1_rt10.2 <- df_pca_rq1(shared_comp_normalized_rt10.2)
# df_pca_rq1_rt10.3 <- df_pca_rq1(shared_comp_normalized_rt10.3)

# Conduct principal component analysis (PCA):
colnames(df_pca_rq1_rt10.1[[2]]) <- c("Fuel type")
# colnames(df_pca_rq1_rt10.2[[2]]) <- c("Fuel type", "Gas Stations")
# colnames(df_pca_rq1_rt10.3[[2]]) <- c("Fuel type", "Gas Stations")

p_rq1_rt10.1 <- pca(mat = df_pca_rq1_rt10.1[[1]], metadata = df_pca_rq1_rt10.1[[2]])
# p_rq1_rt10.2 <- pca(mat = df_pca_rq1_rt10.2[[1]], metadata = df_pca_rq1_rt10.2[[2]])
# p_rq1_rt10.3 <- pca(mat = df_pca_rq1_rt10.3[[1]], metadata = df_pca_rq1_rt10.3[[2]])

# A bi-plot
library(ggalt)
PCAtools::biplot(p_rq1_rt10.1,
                 lab = NULL, # #row.names(p_rq1_rt10.2$metadata)
                 colby = 'Fuel type',
                 hline = 0, vline = 0,
                 legendPosition = 'right', labSize = 5,
                 sizeLoadingsNames = 5,
                 showLoadings = FALSE,
                 # showLoadingsNames = FALSE,
                 ntopLoadings = 10,
                 pointSize = 4, 
                 legendLabSize = 15,
                 legendTitleSize = 16,
                 legendIconSize = 6) + coord_fixed(ratio = 1) 

# Retrieve compound name of top 100 loading
loadingS_rq1_sorted <- p_rq1$loadings %>% arrange(PC1, PC2)
toploadings_rq1 <- rownames(loadingS_rq1_sorted[c(1:50, nrow(loadingS_rq1_sorted):(nrow(loadingS_rq1_sorted) - 50)),1:2])

# Cross check similarity between top 100 loading and RQ1 Wilcoxon test result ASTM_alpha0.1 - 12 compounds
idx <- which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, PASS WILCOXON TEST (data imputed with LOD), alpha threshold < 0.1 (rt1thres = 0.2)`, 
                           "x"))

name <- unlist(flatten(ASTM_list[idx, 16]))

ASTM_rq1_cat5_wilcoxon_stats_alpha0.1 <- gsub("x_", "", name)

intersect_PCA_Wilcoxon_alpha0.1 <- intersect(toploadings_rq1, ASTM_rq1_cat5_wilcoxon_stats_alpha0.1) # 15 ASTM compounds similar between top 200 loading PCA and ASTM_alpha0.1

# Pairs plot
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


# explore further the collapsed_compounds that are driving these differences along each PC.
plotloadings(p_rq1,
             rangeRetain = 0.05, # top 5% variables = top/bottom 5% of the loadings range per PC
             caption = 'Top 10% variables',
             labSize = 4)


eigencorplot(p, metavars = c('fuel_type', 'gas_station'))




# Research Question 2 PCA ##########################
# Category 2: Compounds that have >=2 data points for each gas stations
df_pca_rq2 <- function(data) {
  mydata2 <- data %>%
    filter(., fuel_type %in% "Gas") %>%
    select(sample_name, collapsed_compound, Percent_Area) %>%
    mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
    mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
    # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(sample_name, collapsed_compound) %>%
    summarise(across(Percent_Area, mean)) %>%
    pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
    column_to_rownames(., var = "collapsed_compound")
  
  # transpose the rows and columns
  transpose_mydata2 <- data.table::transpose(mydata2)
  rownames(transpose_mydata2) <- colnames(mydata2)
  colnames(transpose_mydata2) <- rownames(mydata2)
  transpose_mydata2 <- transpose_mydata2 %>%
    rownames_to_column(., var = "sample_name")
  
  transpose_mydata2_new <- transpose_mydata2 %>%
    # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_6", 
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_4", 
                                              ifelse(str_detect(sample_name, "F005"), "Station_3", 
                                                     ifelse(str_detect(sample_name, "F003"), "Station_2", 
                                                            ifelse(str_detect(sample_name, "F008"), "Station_5", "Composite"))))))) %>%
    relocate(gas_station, .after = sample_name)
  
  # Category 1 (rq2_cat1) : Compound found in only 1 gas station and not in any other
  rq2_cat1 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))
  # Category 2 (rq2_cat2): Compound found in >=2 gas stations
  rq2_cat2 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))
  
  # if compounds has only 1 record
  rq2_cat1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 1]
  # if compounds has 2 record
  temp_1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 2]
  rq2_cat2_col_id <- c()
  rq2_cat1_col_id <- c()
  for (col in 1:ncol(temp_1)) {
    idx1 <- which(!is.na(temp_1[,col]))
    # if 2 records from 2 different gas stations
    if (length(unique(transpose_mydata2_new[idx1]$gas_station)) > 1) {
      rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
    } else { # if 2 records from same gas stations
      rq2_cat1_col_id <- c(rq2_cat1_col_id, col)
    }
  }
  rq2_cat1 <- base::cbind(rq2_cat1, temp_1[, rq2_cat1_col_id])
  rq2_cat2 <- temp_1[, rq2_cat2_col_id]
  # if compounds have 3 records
  temp_2 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 3]
  rq2_cat2_col_id <- c()
  rq2_cat1_col_id <- c()
  for (col in 1:ncol(temp_2)) {
    idx2 <- which(!is.na(temp_2[,col]))
    # if 3 records from different gas stations
    if (length(unique(transpose_mydata2_new[idx2,]$gas_station)) > 1) {
      rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
    } else {  # if 3 records from same gas stations
      rq2_cat1_col_id <- c(rq2_cat1_col_id, col)
    }
  }
  rq2_cat1 <- base::cbind(rq2_cat1, temp_2[, rq2_cat1_col_id])
  rq2_cat2 <- base::cbind(rq2_cat2, temp_2[, rq2_cat2_col_id])
  
  # if compounds have >=4 records -> it definitely appear in >=2 Gas stations (712 compounds)
  temp_3 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) > 3]
  rq2_cat2 <- base::cbind(rq2_cat2, temp_3)
  
  # Insert sample name and gas station codes to Cat1 and Cat2 data frames ====================================
  rq2_cat1 <- base::cbind(rq2_cat1, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  rq2_cat2 <- base::cbind(rq2_cat2, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  # Subset compounds that have >= 2 obs per gas station for Wilcoxon test ------------------------------------------------
  temp_4 <- base::cbind(temp_3, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  col_id <- c()
  
  for (col in 3:ncol(temp_4)) {
    record_count <- c()
    idx <- which(!is.na(temp_4[,col]))
    
    for (station in unique(temp_4[idx, "gas_station"])) {
      count <- sum(temp_4[idx,]$gas_station == station)
      record_count <- c(record_count, count) 
    }
    
    if (sum(record_count >= 2) >= 2) {
      col_id <- c(col_id, col)
    } else {
      next
    }
  }
  
  # Data frame of compounds that have >= 2 obs per gas station
  rq2_cat2_stats <- temp_4[, col_id] # 506 compounds
  
  # Impute NA with LOD  for RQ1 Category 2: compounds in >=2 Gas stations
  for (col in 1:ncol(rq2_cat2_stats)) { 
    rq2_cat2_stats[which(is.na(rq2_cat2_stats[,col])), col] <- runif(length(which(is.na(rq2_cat2_stats[,col]))), 
                                                                     min = sort(data$Percent_Area)[1],
                                                                     max = sort(data$Percent_Area)[2])
  }
  
  
  rq2_cat2_stats_imputed <- base::cbind(rq2_cat2_stats, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  # transpose the rows and columns
  rq2_pca <- rq2_cat2_stats_imputed[,-2] %>% column_to_rownames(., var = "sample_name")
  df_X_rq2 <- data.table::transpose(rq2_pca)
  rownames(df_X_rq2) <- colnames(rq2_pca)
  colnames(df_X_rq2) <- rownames(rq2_pca)
  
  # table for information ((rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq2 <- data.frame(colnames(df_X_rq2)) 
  colnames(metadata_X_rq2) <- c('sample_name')
  
  metadata_X_rq2 <- metadata_X_rq2 %>%
    # # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_6", 
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_4", 
                                              ifelse(str_detect(sample_name, "F005"), "Station_3", 
                                                     ifelse(str_detect(sample_name, "F003"), "Station_2", 
                                                            ifelse(str_detect(sample_name, "F008"), "Station_5", "Composite"))))))) %>%
    column_to_rownames(., var = "sample_name") 
  
  # check that sample names match exactly between pdata and expression data 
  print(all(colnames(df_X_rq2) == rownames(metadata_X_rq2)))
  
  return(list(df_X_rq2, metadata_X_rq2))
}  

df_pca_rq2_rt10.1 <- df_pca_rq2(shared_comp_rt10.1)
# df_pca_rq2_rt10.2 <- df_pca_rq2(shared_comp_normalized_rt10.2)
# df_pca_rq2_rt10.3 <- df_pca_rq2(shared_comp_normalized_rt10.3)

# Conduct principal component analysis (PCA):
# Give a better column name for pretty biplot later on
colnames(df_pca_rq2_rt10.1[[2]]) <- c("Gas Stations")
# colnames(df_pca_rq2_rt10.2[[2]]) <- c("Gas Stations")
# colnames(df_pca_rq2_rt10.3[[2]]) <- c("Gas Stations")

p_rq2_rt10.1 <- pca(mat = df_pca_rq2_rt10.1[[1]], metadata = df_pca_rq2_rt10.1[[2]])
# p_rq2_rt10.2 <- pca(mat = df_pca_rq2_rt10.2[[1]], metadata = df_pca_rq2_rt10.2[[2]])
# p_rq2_rt10.3 <- pca(mat = df_pca_rq2_rt10.3[[1]], metadata = df_pca_rq2_rt10.3[[2]])

# A bi-plot
PCAtools::biplot(p_rq2_rt10.1,
                 lab = NULL, #row.names(p_rq2_rt10.1$metadata),
                 colby = 'Gas Stations',
                 hline = 0, vline = 0,
                 legendPosition = 'right',labSize = 5,
                 sizeLoadingsNames = 5,
                 showLoadings = FALSE,
                 # showLoadingsNames = FALSE,
                 ntopLoadings = 8,
                 pointSize = 4, 
                 legendLabSize = 15,
                 legendTitleSize = 16,
                 legendIconSize = 6) + coord_fixed(ratio = 1) 

# Pairs plot
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


# explore further the collapsed_compounds that are driving these differences along each PC
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


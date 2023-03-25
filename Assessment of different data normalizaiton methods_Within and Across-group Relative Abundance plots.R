# Load required packages
library(ggplot2)
library(dplyr)
library(PCAtools)
library(stats)
library(tidyr)

# PCA biplot between gas station =====================
## TSN --------------
pca_rq2_TSN <- function(data) {
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
    # # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_7",
                                              ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                     ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                            ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
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
  metadata_X_rq2 <- data.frame(unique((data %>%
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
  print(all(colnames(df_X_rq2) == rownames(metadata_X_rq2)))
  
  return(list(df_X_rq2, metadata_X_rq2))
}  
df_pca_rq2_rt10.1_TSN <- pca_rq2_TSN(shared_comp_normalized_rt10.1)
colnames(df_pca_rq2_rt10.1_TSN[[2]]) <- c("Gas Stations")
p_rq2_rt10.1_TSN <- pca(mat = df_pca_rq2_rt10.1_TSN[[1]], metadata = df_pca_rq2_rt10.1_TSN[[2]])
biplot(p_rq2_rt10.1_TSN,
       lab = NULL, #row.names(p_rq2_rt10.1$metadata),
       colby = 'Gas Stations',
       hline = 0, vline = 0,
       legendPosition = 'right',labSize = 5,
       sizeLoadingsNames = 5,
       showLoadings = TRUE,
       # showLoadingsNames = FALSE,
       ntopLoadings = 8,
       pointSize = 4, 
       legendLabSize = 15,
       legendTitleSize = 16,
       legendIconSize = 6)

## Log10 --------------
pca_rq2_log10 <- function(data) {
  mydata2 <- data %>%
    filter(., fuel_type %in% "Gas") %>%
    select(sample_name, collapsed_compound, Log_Area) %>%
    mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
    mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
    # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(sample_name, collapsed_compound) %>%
    summarise(across(Log_Area, mean)) %>%
    pivot_wider(names_from = sample_name, values_from = Log_Area) %>%
    column_to_rownames(., var = "collapsed_compound")
  
  # transpose the rows and columns
  transpose_mydata2 <- data.table::transpose(mydata2)
  rownames(transpose_mydata2) <- colnames(mydata2)
  colnames(transpose_mydata2) <- rownames(mydata2)
  transpose_mydata2 <- transpose_mydata2 %>%
    rownames_to_column(., var = "sample_name")
  
  transpose_mydata2_new <- transpose_mydata2 %>%
    # # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_7",
                                              ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                     ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                            ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
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
                                                                     min = sort(data$Log_Area)[1],
                                                                     max = sort(data$Log_Area)[2])
  }
  
  
  rq2_cat2_stats_imputed <- base::cbind(rq2_cat2_stats, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  # transpose the rows and columns
  rq2_pca <- rq2_cat2_stats_imputed[,-2] %>% column_to_rownames(., var = "sample_name")
  df_X_rq2 <- data.table::transpose(rq2_pca)
  rownames(df_X_rq2) <- colnames(rq2_pca)
  colnames(df_X_rq2) <- rownames(rq2_pca)
  
  # table for information ((rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq2 <- data.frame(unique((data %>%
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
  print(all(colnames(df_X_rq2) == rownames(metadata_X_rq2)))
  
  return(list(df_X_rq2, metadata_X_rq2))
}  
df_pca_rq2_rt10.1_log10 <- pca_rq2_log10(shared_comp_normalized_rt10.1)
colnames(df_pca_rq2_rt10.1_log10[[2]]) <- c("Gas Stations")
p_rq2_rt10.1_log10 <- pca(mat = df_pca_rq2_rt10.1_log10[[1]], metadata = df_pca_rq2_rt10.1_log10[[2]])
biplot(p_rq2_rt10.1_log10,
       lab = NULL, #row.names(p_rq2_rt10.1$metadata),
       colby = 'Gas Stations',
       hline = 0, vline = 0,
       legendPosition = 'right',labSize = 5,
       sizeLoadingsNames = 5,
       showLoadings = TRUE,
       # showLoadingsNames = FALSE,
       ntopLoadings = 8,
       pointSize = 4, 
       legendLabSize = 15,
       legendTitleSize = 16,
       legendIconSize = 6)

## boxcox --------------
pca_rq2_boxcox <- function(data) {
  mydata2 <- data %>%
    filter(., fuel_type %in% "Gas") %>%
    select(sample_name, collapsed_compound, boxcox_area) %>%
    mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
    mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
    # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
    group_by(sample_name, collapsed_compound) %>%
    summarise(across(boxcox_area, mean)) %>%
    pivot_wider(names_from = sample_name, values_from = boxcox_area) %>%
    column_to_rownames(., var = "collapsed_compound")
  
  # transpose the rows and columns
  transpose_mydata2 <- data.table::transpose(mydata2)
  rownames(transpose_mydata2) <- colnames(mydata2)
  colnames(transpose_mydata2) <- rownames(mydata2)
  transpose_mydata2 <- transpose_mydata2 %>%
    rownames_to_column(., var = "sample_name")
  
  transpose_mydata2_new <- transpose_mydata2 %>%
    # # Grouping samples into respective Gas stations
    mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                                ifelse(str_detect(sample_name, "F001"), "Station_1",
                                       ifelse(str_detect(sample_name, "F007"), "Station_7",
                                              ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                     ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                            ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
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
                                                                     min = sort(data$boxcox_area)[1],
                                                                     max = sort(data$boxcox_area)[2])
  }
  
  
  rq2_cat2_stats_imputed <- base::cbind(rq2_cat2_stats, transpose_mydata2_new[, 1:2]) %>%
    relocate(sample_name, gas_station, .before = everything())
  
  # transpose the rows and columns
  rq2_pca <- rq2_cat2_stats_imputed[,-2] %>% column_to_rownames(., var = "sample_name")
  df_X_rq2 <- data.table::transpose(rq2_pca)
  rownames(df_X_rq2) <- colnames(rq2_pca)
  colnames(df_X_rq2) <- rownames(rq2_pca)
  
  # table for information ((rows are sample IDs, columns are sample information) -----------------------
  metadata_X_rq2 <- data.frame(unique((data %>%
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
  print(all(colnames(df_X_rq2) == rownames(metadata_X_rq2)))
  
  return(list(df_X_rq2, metadata_X_rq2))
}  
df_pca_rq2_rt10.1_boxcox <- pca_rq2_boxcox(shared_comp_normalized_rt10.1)
colnames(df_pca_rq2_rt10.1_boxcox[[2]]) <- c("Gas Stations")
p_rq2_rt10.1_boxcox <- pca(mat = df_pca_rq2_rt10.1_boxcox[[1]], metadata = df_pca_rq2_rt10.1_boxcox[[2]])
biplot(p_rq2_rt10.1_boxcox,
       lab = NULL, #row.names(p_rq2_rt10.1$metadata),
       colby = 'Gas Stations',
       hline = 0, vline = 0,
       legendPosition = 'right',labSize = 5,
       sizeLoadingsNames = 5,
       showLoadings = TRUE,
       # showLoadingsNames = FALSE,
       ntopLoadings = 8,
       pointSize = 4, 
       legendLabSize = 15,
       legendTitleSize = 16,
       legendIconSize = 6)



# HCA biplot ==========================

## TSN
stats::hclust()
## Log10
## boxcox


# Calculate relative abundances within each group ==============
abundance_data_rel <- abundance_data %>%
  group_by(Group) %>%
  mutate(rel_abundance_original = Area / sum(Area))

# Plot within-group relative abundances
ggplot(abundance_data_rel, aes(x = Group, y = rel_abundance)) +
  geom_boxplot() +
  ggtitle("Within-Group Relative Abundances") +
  ylab("Relative Abundance")

# Calculate relative abundances across all groups ==============
abundance_data_rel_total <- abundance_data_rel %>%
  group_by(Sample) %>%
  mutate(rel_abundance_total = Abundance / sum(Abundance))

# Plot across-group relative abundances
ggplot(abundance_data_rel_total, aes(x = Sample, y = rel_abundance_total)) +
  geom_boxplot() +
  ggtitle("Across-Group Relative Abundances") +
  ylab("Relative Abundance")


# P-value histograms =============
## TSN
## Log10
## boxcox

# Histogram of original data, TSN, log10 and Boxcox normalized data -----------------------

# original data
ggplot(data = shared_comp_normalized_rt10.1) +
  geom_histogram(aes(x=Area), 
                 bins = 2000000/20000) +
  scale_x_continuous(limits = c(0, 60000000)) +
  scale_y_continuous(limits = c(0, 6000)) +
  labs(x = "original Peak area") +
  theme_classic(base_size = 20)

# TSN data
ggplot(data = shared_comp_normalized_rt10.1) +
  geom_histogram(aes(x=Percent_Area),
                 bins = 100) +
  scale_x_continuous(limits = c(0, 0.04)) +
  scale_y_continuous(limits = c(0, 6000)) +
  labs(x = "Total Sum normalized area") +
  theme_classic(base_size = 20)

# log10 data
ggplot(data = shared_comp_normalized_rt10.1) +
  geom_histogram(aes(x=Log_Area),
                 bins = 100) +
  labs(x = "Log10 normalized area") +
  theme_classic(base_size = 20)

# Boxcox data
ggplot(data = shared_comp_normalized_rt10.1) +
  geom_histogram(aes(x=boxcox_area), 
                 bins = 100) +
  labs(x = "Box-Cox normalized area") +
  theme_classic(base_size = 20)



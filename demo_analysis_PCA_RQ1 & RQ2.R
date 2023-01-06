# Loading Packages --------------------------------------------------------
# library(DiagrammeR)
library(ggplot2)
library(readxl)
library(ggpubr)
library(gtable)
library(gridExtra)
library(viridis)
library(wesanderson)
library(tidyverse)
library(lubridate)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(data.table)
library(corrplot)
library(purrr)
library(stringr)
library(stringi)
library(grid)
library(missMDA)
library(tsne)
library(Rtsne)
library(rgl)
library(collapse)
library(plotly)
library(umap)
library(sqldf)
library(multiway)
library(mdatools)
library(writexl)
library(ggsignif)

# Functions -------------------------------------------------------------------------------------------------------
# Filtering matched compound names
filtering <- function(df, filter_list) {
  clean_data <-  copy(df)
  for (ele in filter_list) {
    clean_data <- clean_data %>%
      filter(!grepl(ele, Compound))
  }
  return(clean_data)
}

# Filtering limit of observations
limit_obser <- function(df_list, file_list, cap) {
  df_list_filter_area <- list()
  df_list_removed_area <- list()
  for (i in 1:length(df_list)) {
    df_list_filter_area[[i]] <- df_list[[i]] %>%
      filter(., Area > cap) %>%
      mutate(sample_name = file_list[[i]]) %>%
      # Grouping samples into fuel types
      mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                       ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))) %>%
      # Grouping samples into respective Gas stations
      mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9", 
                                  ifelse(str_detect(sample_name, "F001"), "Station_1",
                                         ifelse(str_detect(sample_name, "F007"), "Station_7", 
                                                ifelse(str_detect(sample_name, "F005"), "Station_5", 
                                                       ifelse(str_detect(sample_name, "F003"), "Station_3", 
                                                              ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite")))))))
    
    df_list_removed_area[[i]] <- df_list[[i]] %>%
      filter(., Area <= cap) %>%
      mutate(sample_name = file_list[[i]]) %>%
      # Grouping samples into fuel types
      mutate(fuel_type = ifelse(str_detect(sample_name, "DieselComp"), "DieselComp", 
                                ifelse(str_detect(sample_name, "GasComp"), "GasComp",
                                       ifelse(str_detect(sample_name, "D"), "Diesel", "Gas")))) %>%
      # Grouping samples into respective Gas stations
      mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9", 
                                  ifelse(str_detect(sample_name, "F001"), "Station_1",
                                         ifelse(str_detect(sample_name, "F007"), "Station_7", 
                                                ifelse(str_detect(sample_name, "F005"), "Station_5", 
                                                       ifelse(str_detect(sample_name, "F003"), "Station_3", 
                                                              ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite")))))))
  }
  return(list(df_list_filter_area, df_list_removed_area))
}

# Notin function
`%notin%` <- Negate(`%in%`)

# Grouping compounds based on RT1, RT2, and Ion1 - Version 1
grouping_comp_ver1 <- function(data, rt1thres, rt2thres, ion1thres, ion2thres) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rt1thres: RT1 threshold window
  # rt2thres: RT2 threshold window
  # ion1thres: Ion1 threshold window
  # ion2thres: Ion2 threshold window
  # region_applied: list of x-y coordinate for different regions to applied different threshold window, x-axis is RT1, y-axis is RT2
  dat <- copy(data)
  
  # Initialize the compound column filled with NA values
  dat$collapsed_compound <- NA
  i <- 1
  
  for (row in 1:nrow(dat)) {
    # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    rt1 <- dat[row,]$RT1
    rt2 <- dat[row,]$RT2
    ion1 <- dat[row,]$Ion1
    ion2 <- dat[row,]$Ion2
    
    idx_thres <- which(dat$RT1 <= (rt1 + rt1thres) & dat$RT1 >= (rt1 - rt1thres) & 
                         dat$RT2 <= (rt2 + rt2thres) & dat$RT2 >= (rt2 - rt2thres) & 
                         dat$Ion1 <= (ion1 + ion1thres) & dat$Ion1 >= (ion1 - ion1thres) & 
                         dat$Ion2 <= (ion2 + ion2thres) & dat$Ion2 >= (ion2 - ion2thres) &
                         is.na(dat$collapsed_compound))
    
    if (identical(idx_thres, integer(0))) {
      next
    }
    else {
      dat[idx_thres, "collapsed_compound"] <- paste0("Compound_", i, ".")
      i <- i + 1
    }  
  }
  
  return(dat)
}

# Grouping compounds based on RT1, RT2, and Ion1 - Version 2
grouping_comp_ver2 <- function(data) {
  
  # create empty list, each sub-list is a compound group with following criteria:
  # rt1thres: RT1 threshold window
  # rt2thres: RT2 threshold window
  # ion1thres: Ion1 threshold window
  # ion2thres: Ion2 threshold window
  # region_applied: list of x-y coordinate for different regions to applied different threshold window, x-axis is RT1, y-axis is RT2
  data <- copy(all_data_pre_norm_filter_area)
  
  region_num <- as.numeric(base::readline("Please input the number of dividing region: "))
  rt1thres <- as.numeric(base::readline("Please input the RT1 window threshold for the region applied: "))
  rt2thres <- as.numeric(base::readline("Please input the RT2 window threshold for the region applied: "))
  ion1thres <- as.numeric(base::readline("Please input the Ion1 window threshold for the region applied: "))
  ion2thres <- as.numeric(base::readline("Please input the Ion2 window threshold for the region applied: "))
  
  # Initialize the compound column filled with NA values
  data$compound <- NA
  i <- 1
  
  for (reg_num in 1:region_num) {
    # User will input the coordinate of region that they want to applied a specific threshold of RT1
    region_x1 <- as.numeric(base::readline("Please input the bottom-left coordinate of the region applied: "))
    region_x2 <- as.numeric(base::readline("Please input the bottom-right coordinate of the region applied: "))
    region_y1 <- as.numeric(base::readline("Please input the top-left coordinate of the region applied: "))
    region_y2 <- as.numeric(base::readline("Please input the top-right coordinate of the region applied: "))
    
    region_applied <- c(region_x1, region_x2, region_y1, region_y2)
    
    # all_data_pre_norm_filter_area <- bind_rows(list_remaining_area) %>%
    #   arrange(RT1, RT2)
    # 
    # rt1thres <- 0.2
    # rt2thres <- 0.2
    # ion1thres <- 0.05
    # region_applied <- c(2, 72, 1, 6)
    
    # Filter data frame so that it only contain data in the region applied 
    idx_region <- which(data$RT1 >= region_applied[1] & data$RT1 <= region_applied[2] & 
                          data$RT2 >= region_applied[3] & data$RT2 <= region_applied[4])
    
    for (row in idx_region) {
      # filter data by index, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
      rt1 <- data[idx_region,][row,]$RT1
      rt2 <- data[idx_region,][row,]$RT2
      ion1 <- data[idx_region,][row,]$Ion1
      ion2 <- data[idx_region,][row,]$Ion2
      
      idx_thres <- which(data[idx_region,]$RT1 <= (rt1 + rt1thres) & data[idx_region,]$RT1 >= (rt1 - rt1thres) & 
                           data[idx_region,]$RT2 <= (rt2 + rt2thres) & data[idx_region,]$RT2 >= (rt2 - rt2thres) & 
                           data[idx_region,]$Ion1 <= (ion1 + ion1thres) & data[idx_region,]$Ion1 >= (ion1 - ion1thres) & 
                           data[idx_region,]$Ion2 <= (ion2 + ion2thres) & data[idx_region,]$Ion2 >= (ion2 - ion2thres) &
                           is.na(data[idx_region,]$compound))
      
      if (identical(idx_thres, integer(0))) {
        next
      }
      else {
        data[idx_region,][idx_thres, "compound"] <- paste0("Compound_", i, ".")
        i <- i + 1
      }  
    }
  }
  
  return(data)
}

# Filtering similar and unique compound
comp_filter_ver1 <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$collapsed_compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$sample_name)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$sample_name)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# Filtering similar and unique compound
comp_filter_ver2 <- function(data, n) {
  all_similar_compounds_idx <- c()
  all_other_compounds_idx <- c()
  all_unique_compounds_idx <- c()
  
  for (comp_grp in unique(data$collapsed_compound)) {
    # filter data by indexing, ALWAYS DO THIS INSTEAD OF CREATE SUBSET DATAFRAME
    idx <- which(grepl(comp_grp, data$collapsed_compound, fixed = TRUE))
    
    if (length(unique(data[idx,]$gas_station)) > (n - 1)) {
      all_similar_compounds_idx <- c(all_similar_compounds_idx, idx)
    }
    else if (length(unique(data[idx,]$gas_station)) < 2) {
      all_unique_compounds_idx <- c(all_unique_compounds_idx, idx)
    }
    else {
      all_other_compounds_idx <- c(all_other_compounds_idx, idx)
    }
  }
  return(list(all_similar_compounds_idx, all_other_compounds_idx, all_unique_compounds_idx))
}

# Probabilistic Quotient Normalization
pqn <- function(X, n = "median", QC = NULL) {
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  if (!is.null(QC)) {
    # if QC vector exists, use this as reference spectrum
    if (length(QC) == 1) {
      # only 1 reference sample given
      mX <- as.numeric(X[QC, ])
    } else {
      if (n == "mean") {
        mX <- as.numeric(colMeans(X[QC, ]))
      }
      if (n == "median") {
        mX <- as.numeric(apply(X[QC, ], 2, median))
      }
    }
  } else {
    # otherwise use the mean or median of all samples as reference sample
    if (n == "mean") {
      mX <- as.numeric(colMeans(X))
    }
    if (n == "median") {
      mX <- as.numeric(apply(X, 2, median))
    }
  }
  
  # do the actual normalisation
  for (a in 1:nrow(X)) {
    X.norm[a, ] <- as.numeric(X[a, ] / median(as.numeric(X[a, ] / mX)))
  }
  
  return(X.norm)
}

# Relative log abundance plots
RlaPlots <- function(inputdata, type=c("ag", "wg"), cols=NULL,
                     cex.axis=0.8, las=2, ylim=c(-2, 2), oma=c(7, 4, 4, 2) + 0.1, ...) {
  type <- match.arg(type)
  groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
  unique.groups <- levels(groups)
  if (is.null(cols)) 
    cols <- ColList(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(inputdata))))
  for (ii in 1:length(inputdata[, 1])) 
    box_cols[ii] <- cols[which(unique.groups == inputdata[, 1][ii])]
  
  # Within groups
  if(type == "wg") {
    out_data<-data.frame()
    for (grp in unique.groups) {
      submat <- inputdata[which(inputdata[, 1] == grp), -1]
      med_vals <- apply(submat, 2, median)
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_data <- rbind(out_data, swept_mat)
    }
    # Across groups (i.e. type == "ag")
  } else  {
    med_vals <- apply(inputdata[, -1], 2, median)
    out_data <- sweep(inputdata[, -1], 2, med_vals, "-")
  }
  
  boxplot(t(out_data),
          cex.axis=cex.axis,                 # font size
          las=las,                           # label orientation
          col=box_cols,                      # colours
          ylim=ylim,                         # y-axis range
          oma=oma,                           # outer margin size
          ...
  )
  
  abline(h=0)
}


# STEP 1.1: Data import --------------------------------------------
file_list <- list.files(pattern = '*.xlsx')

# Pipe operator for isolating IL types
indi_IL_file_list <- file_list %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")] %>%
  .[!str_detect(., "Station")] %>%
  .[!str_detect(., "Workflow")] %>%
  .[!str_detect(., "ASTM")] %>%
  .[!str_detect(., "data_table")]

gas_only_file_list <- file_list %>%
  .[!str_detect(., "D")]  %>%
  .[!str_detect(., "DieselComp")] %>%
  .[!str_detect(., "GasComp")] %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")] %>%
  .[!str_detect(., "Station")] %>%
  .[!str_detect(., "Workflow")] %>%
  .[!str_detect(., "ASTM")] %>%
  .[!str_detect(., "data_table")]

non_gas_file_list <- file_list[!file_list %in% gas_only_file_list] %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")] %>%
  .[!str_detect(., "Station")] %>%
  .[!str_detect(., "Workflow")] %>%
  .[!str_detect(., "ASTM")] %>%
  .[!str_detect(., "data_table")]

gas_diesel_only_file_list <- file_list %>%
  .[!str_detect(., "DieselComp")] %>%
  .[!str_detect(., "GasComp")] %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")] %>%
  .[!str_detect(., "Station")] %>%
  .[!str_detect(., "Workflow")] %>%
  .[!str_detect(., "ASTM")] %>%
  .[!str_detect(., "data_table")]


# Import IL samples to list
df_list <- purrr::map(indi_IL_file_list, read_xlsx, sheet = "Results")

df <- dplyr::bind_rows(df_list)

# remove spaces in column names 
for (i in 1:length(df_list)) {
  colnames(df_list[[i]]) <- gsub(" ", "", colnames(df_list[[i]]))
}

# STEP 1.2A: Filtering out column bleed, solvent and BTEX and minimum area observations--------------------------------------
df_list_clean <- purrr::map(df_list, filtering, filter_list = c("^Carbon disulfide$", 
                                                                "Cyclotrisiloxane..hexamethyl",
                                                                "Cyclotetrasiloxane..octamethyl",
                                                                "^Benzene$",
                                                                "^Toluene$",
                                                                "^Ethylbenzene$",
                                                                "Xylene")) 
df_clean <- dplyr::bind_rows(df_list_clean)
# STEP 1.2B Filtering out limit of observations----------------------------
list_remaining_area <- limit_obser(df_list_clean, indi_IL_file_list, cap = 50000)[[1]]
list_removed_area <- limit_obser(df_list_clean, indi_IL_file_list, cap = 50000)[[2]]
# STEP 1.3: Grouping compounds based on RT1, RT2, Ion1 -----------------------------------------------------------------------
# STEP 1.3A: Generate 1 grand data frame of all 31 IL samples

all_data_pre_norm_filter_area <- bind_rows(list_remaining_area) %>%
  arrange(RT1, RT2)

# STEP 1.3B: Grouping compounds based on RT1, RT2, Ion1 threshold ----------------------------------------
all_data_pre_norm_grouped_filter_area <- grouping_comp_ver1(all_data_pre_norm_filter_area,
                                                            rt1thres = 0.2, # rt1thres = 0.2
                                                            rt2thres = 0.125,
                                                            ion1thres = 0.05, # Ion 1 and 2 indicates molecular structure (2 most prevalent mass-to-charge)
                                                            ion2thres = 0.05) 

# STEP 2: Identify shared and unique compound groups across samples ------------------------------------------------
idx_list_filter_area_samples <- comp_filter_ver1(all_data_pre_norm_grouped_filter_area, 
                                                 length(indi_IL_file_list))

similar_compounds_filter_area_samples <- all_data_pre_norm_grouped_filter_area[idx_list_filter_area_samples[[1]],] 
other_compounds_filter_area_samples <- all_data_pre_norm_grouped_filter_area[idx_list_filter_area_samples[[2]],] 
unique_compounds_filter_area_samples <- all_data_pre_norm_grouped_filter_area[idx_list_filter_area_samples[[3]],]

# Combine similar_compounds_filter_area and other_compounds_filter_area to one data frame 
shared_comp_samples <- bind_rows(similar_compounds_filter_area_samples, other_compounds_filter_area_samples)

# STEP 3: Data Normalization ================================================================================
add_data_normalization <- function(data) {
  temp_list <- list()
  i <- 1
  for (sample in unique(data$sample_name)) {
    df <- data[which(data$sample_name == sample),] %>%
      # Log-based normalization
      mutate(Log_Area = log10(Area)) %>%
      mutate(Log_Height = log10(Height)) %>%
      # TSN - Percent-based normalization
      mutate(Percent_Area = Area/sum(.$Area)) %>%
      mutate(Percent_Height = Height/sum(.$Height))
    temp_list[[i]] <- df
    i <- i + 1
  }
  # Then combine data again to 1 grand data frame
  newdata <- dplyr::bind_rows(temp_list)
  return(newdata)
}

shared_comp_normalized <- add_data_normalization(shared_comp_samples)

# STEP 4: Visualization with PCA(-HCPC)

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

cat_5 <- df_X_rq1[(rowSums(!is.na(df_X_rq1[, 1:21])) >= 2) & 
                       (rowSums(!is.na(df_X_rq1[, 22:25])) >= 2),] # -> 374 out of 4769 compounds

for (r in 1:nrow(cat_5)) { 
  cat_5[r, which(is.na(cat_5[r,]))] <- runif(length(which(is.na(cat_5[r,]))),
                                             min = sort(shared_comp_normalized$Percent_Area)[1],
                                             max = sort(shared_comp_normalized$Percent_Area)[2])
}

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

all(colnames(cat_5) == rownames(metadata_X_rq1))

# Conduct principal component analysis (PCA): ---------------------------------
library(PCAtools)

p_rq1 <- pca(cat_5, metadata = metadata_X_rq1) # Classify on sample

# A bi-plot -------------
biplot(p_rq1,
       lab = row.names(p_rq1$metadata),
       colby = 'fuel_type',
       hline = 0, vline = 0,
       legendPosition = 'right', labSize = 4,
       sizeLoadingsNames = 4,
       showLoadings = TRUE,
       ntopLoadings = 50)

# Retrieve compound name of top 100 loading
loadingS_rq1_sorted <- p_rq1$loadings %>% arrange(PC1, PC2)
toploadings_rq1 <- rownames(loadingS_rq1_sorted[c(1:100, nrow(loadingS_rq1_sorted):(nrow(loadingS_rq1_sorted) - 100)),1:2])

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
rq2_pca <- rq2_cat2_stats[,-2] %>% column_to_rownames(., var = "sample_name")
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
library(PCAtools)

p_rq2 <- pca(df_X_rq2, metadata = metadata_X_rq2)

# A bi-plot -------------
biplot(p_rq2,
       lab = row.names(p_rq2$metadata),
       colby = 'gas_station',
       hline = 0, vline = 0,
       legendPosition = 'right', labSize = 5,
       sizeLoadingsNames = 3,
       showLoadings = TRUE,
       ntopLoadings = 50)

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


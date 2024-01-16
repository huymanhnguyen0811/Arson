# ASTM compound list
ASTM_list <- read_xlsx(paste0(getwd(), "/ASTM Compound List.xlsx"), sheet = "Sheet1")

# Step 1.1 Data Import
sum(ASTM_list$`step 1.1` == "x") # 71/80 ASTM compounds found

# Step 1.2 Filtering column bleed / solvent / BTEX and minimum area & 1.3 Collapsing compounds based on RT1 RT2 Ion1 Ion2
sum(ASTM_list$`step 1.2 + step 1.3` == "x") # 66/80 ASTM compounds found; because remove 5 BTEX compounds

# Step 2 Sub-setting non-unique compounds to new data frame and Step 3 Data normalization
sum(str_detect(ASTM_list$`step 2 + step 3 (Shared_comp) (rt1_thres = 0.2)`, "x")) # 65/80 ASTM compounds found; lost 1,4-Diethylbenzene

# Step 4A. Research Question 1 ============================================================================================================================
mydata_ASTM <- shared_comp_normalized %>%
  filter(., fuel_type %in% c("Gas", "Diesel")) %>%
  select(sample_name, collapsed_compound, Percent_Area) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(sample_name, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
  relocate(`0220F001D.xlsx`, `0220F009D.xlsx`, `0220F009-2D.xlsx`, `0220F005D.xlsx`, .after = last_col()) %>%
  column_to_rownames(.,var = "collapsed_compound")

ASTM_list[which(str_detect(ASTM_list$`RQ1: Compounds found in only Gas and Diesel`, "x")), "Compound"]

# category 1: Compound with No Diesel records, how we determine the significant compounds ------------------------------
# How many compounds has no diesel records? -> 1137 out of 4769 compounds 
cat_1 <- mydata_ASTM[rowSums(is.na(mydata_ASTM[, c(22:25)])) == 4,] # BEWARE: these compounds may only appear in 2 gas samples

dim(cat_1)

# How many of them are ASTM compounds? - 16 compounds
sum(str_detect(ASTM_list$`RQ1: Category 1: compounds only in Gas`, "x"))

# How many observations - 5284 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% rownames(cat_1)))


ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 1: compounds only in Gas`, "x")), "Compound"]

# Octane (n-alkanes)
# 3-Ethyltoluene (Castle Group)
# 4-Ethyltoluene (Castle Group)
# 1,2,4-Trimethylbenzene
# 1-Ethyl-3,5-Dimethylbenzene (Gang of Four)
# Undecane (n-alkanes)
# Dodecane (n-alkanes)
# Hexylbenzene
# Tridecane (n-alkanes)
# 2,3-Dimethylnaphthalene (Five Fingers)
# 2,6- & 1,3- & 1,7-Dimethylnaphthalene (Five Fingers)
# 1,2,3-Trimethylbenzene
# 1-Ethyl-2,5-dimethylbenzene (2-1)
# 1-Ethyl-2,6-dimethylbenzene 
# 1,2,3,5-Tetramethylbenzene (Tetris)
# 2-Methylindane


# category 2: Compound with No Gas records, how we determine the significant compounds ---------------------------------
# How many compounds has no diesel records? -> 2812 out of 4769 compounds 
cat_2 <- mydata_ASTM[rowSums(is.na(mydata_ASTM[, c(1:21)])) == 21,] 

# How many of them are ASTM compounds? - 14 ASTM compounds
sum(str_detect(ASTM_list$`RQ1: Category 2: compounds only in Diesel`, "x"))

# How many observations - 10225 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% rownames(cat_2)))

# What are the ASTM? 
ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 2: compounds only in Diesel`, "x")), "Compound"]

# Undecane (n-alkanes)
# Dodecane (n-alkanes)
# Tridecane (n-alkanes)
# Pentadecane (n-alkanes)
# Heptadecane (n-alkanes)
# Hexadecane (n-alkanes)
# Nonadecane (n-alkanes)
# Octadecane n-alkanes)
# Eicosane (n-alkanes)
# Heneicosane (n-alkanes)
# Tricosane (n-alkanes)
# 1,2,-Diethylbenzene
# 1,2-Dimethylnaphthalene (Five Fingers)
# 1,3,5-Triethylbenzene


# category 3: Compound with >= 1 Gas record and only 1 diesel record, how we determine the significant compounds ---------------------------
# How many compounds has no diesel records? -> 345 out of 4769 compounds 
cat_3 <- mydata_ASTM[(rowSums(is.na(mydata_ASTM[, 22:25])) == 3) &  # only 1 Diesel record
                       (rowSums(!is.na(mydata_ASTM[, 1:21])) >= 1),] # only 1 Gas record

# How many observations - 10225 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% rownames(cat_3)))

# What are the ASTM? # How many of them are ASTM compounds? - 6 compounds
ASTM_list[which(str_detect(ASTM_list$`RQ1 Category 3: Compounds with only 1 Diesel and at least 1 Gas record`, "x")), "Compound"]

# Decane
# 4-Isopropyltoluene
# 1-Ethyl-3,5-Dimethylbenzene
# Dodecane
# 2,3-Dimethylnaphthalene
# 2-Methylindane


# category 4: only 1 Gas and >=2 Diesel --------------------------------------------------------
cat_4 <- mydata_ASTM[(rowSums(!is.na(mydata_ASTM[, 1:21])) == 1) & 
                       (rowSums(!is.na(mydata_ASTM[, 22:25])) >= 2),] # -> 101 out of 4769 compounds

# How many observations - 634 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% rownames(cat_4)))

# How many of them are ASTM compounds? - 1 compounds
ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 4: Compounds with only one Gas and  >=2 diesel records`, "x")), "Compound"]
# 4-Isopropyltoluene

# category 5: >=2 gas and >=2 diesel (multiple testing) --------------------------------------------------------
cat_5 <- mydata_ASTM[(rowSums(!is.na(mydata_ASTM[, 1:21])) >= 2) & 
                       (rowSums(!is.na(mydata_ASTM[, 22:25])) >= 2),] # -> 374 out of 4769 compounds


# How many observations - 6389 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% rownames(cat_5)))

# How many of them are ASTM compounds?
# BEFORE WILCOXON TEST - 44 ASTM compounds
cat5_b4wilcoxon <- c(ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, BEFORE WILCOXON TEST (data imputed with LOD)`, 
                                                "x")), "Compound"])[[1]]

# cat5_b4wilcoxon <- c("1,8-Dimethylnaphthalene (Five fingers)", "Dimethylindane (24.1321)" , "Dimethylindane (24)", "Dimethylindane (23.8366)", "Pentylbenzene" ,
#   "1,2,3,4-Tetramethylbenzene", "5-Methylindane", "1-Methyl-2-tert-butylbenzene", "1,2,4,5-Tetramethylbenzene", "1-Ethyl-2,4-dimethylbenzene (2-1)" ,
#   "1-Ethyl-2,5-dimethylbenzene" ,"2-Propyltoluene", "1,2,3-Trimethylbenzene", "Isopropylbenzene", "2,6- & 1,3- & 1,7-Dimethylnaphthalene (Five Fingers)", 
#   "1,4-Dimethylnaphthalene (Five Fingers)", "1,6-Dimethylnaphthalene (Five Fingers)", "1-Ethylnaphthalene (Five Fingers)", 
#   "2-Ethylnaphthalene (Five Fingers)", "Tetradecane (n-alkanes)", "1-Methylnaphthalene", "2-Methylnaphthalene" ,"Tridecane", 
#   "Naphthalene (PAHs)", "4,7-Dimethylindane", "Isopentylbenzene", "1-Ethyl-2,3-dimethylbenzene", "Undecane", 
#   "3- & 4-Propyltoluene (Gang of Four)", "1,3-Diethylbenzene (Gang of Four)", "Indane", "2-Isopropyltoluene", 
#   "4-Isopropyltoluene", "3-Isopropyltoluene", "Sec-butylbenzene", "Isobutylbenzene", "1,2,4-Trimethylbenzene", 
#   "2-Ethyltoluene (Castle Group)", "1,3,5-Trimethylbenzene", "4-Ethyltoluene", "3-Ethyltoluene", "Propylbenzene (Castle Group)", 
#   "Nonane (n-alkanes)", "Octane (n-alkanes)")

# AFTER WILCOXON TEST, PASS alpha threshold < 0.1 - 21 ASTM compounds
cat5_passwilcoxon_alpha0.1 <- c(ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, PASS WILCOXON TEST (data imputed with LOD), alpha threshold < 0.1 (rt1thres = 0.2)`, 
                                                "x")), "Compound"])[[1]]

# ASTM compounds AFTER WILCOXON TEST, FAIL alpha threshold < 0.1
fail_wilcoxon_alpha0.1 <- setdiff(cat5_b4wilcoxon, cat5_passwilcoxon_alpha0.1)
fail_wilcoxon_alpha0.1

# INtersect between PCA top loadings and compounds that passed Wilcoxon alpha = 0.1 ====================================
ASTM_intersect_rq1_pca_wilcox <- c()
for (ele in intersect_PCA_Wilcoxon_alpha0.1) {
  ASTM_intersect_rq1_pca_wilcox <- c(ASTM_intersect_rq1_pca_wilcox, 
                                     as.character(ASTM_list[which(str_detect(ASTM_list$`RQ1: Category 5: Compounds with >=2 Gas record and >=2 Diesel record, PASS WILCOXON TEST (data imputed with LOD), alpha threshold < 0.1`, 
                                                                             ele)), "Compound"]))
}


# Research Question 2 ============================================================================================================================
mydata2 <- shared_comp_normalized %>%
  filter(., fuel_type %in% "Gas") %>%
  select(sample_name, collapsed_compound, Percent_Area) %>%
  mutate(sample_name = factor(sample_name, levels = c(unique(sample_name)))) %>%
  mutate(collapsed_compound = factor(collapsed_compound, levels = c(unique(collapsed_compound)))) %>%
  # since we have duplicates with different values of the same compound in some samples, we summarize these values by taking the mean of them
  group_by(sample_name, collapsed_compound) %>%
  summarise(across(Percent_Area, mean)) %>%
  pivot_wider(names_from = sample_name, values_from = Percent_Area) %>%
  column_to_rownames(., var = "collapsed_compound")

ASTM_list[which(str_detect(ASTM_list$`RQ2: compounds found only in Gas`, "x")), "Compound"]

# transpose the rows and columns
transpose_mydata2 <- data.table::transpose(mydata2)
rownames(transpose_mydata2) <- colnames(mydata2)
colnames(transpose_mydata2) <- rownames(mydata2)
transpose_mydata2 <- transpose_mydata2 %>%
  rownames_to_column(., var = "sample_name")

# Grouping samples into respective Gas stations
transpose_mydata2_new <- transpose_mydata2 %>%
  mutate(gas_station = ifelse(str_detect(sample_name, "F009"), "Station_9",
                              ifelse(str_detect(sample_name, "F001"), "Station_1",
                                     ifelse(str_detect(sample_name, "F007"), "Station_7",
                                            ifelse(str_detect(sample_name, "F005"), "Station_5",
                                                   ifelse(str_detect(sample_name, "F003"), "Station_3",
                                                          ifelse(str_detect(sample_name, "F008"), "Station_8", "Composite"))))))) %>%
  relocate(gas_station, .after = sample_name)

# Initiate empty data frame for Category 1 (rq2_cat1) : Compound found in only 1 gas station and not in any other
rq2_cat1 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))

# Initiate empty data frame for Category 2 (rq2_cat2): Compound found in >=2 gas stations
rq2_cat2 <- data.frame(matrix(nrow = nrow(transpose_mydata2_new)))

# if compounds has only 1 record
rq2_cat1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 1]
# if compounds has 2 record
temp_1 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) == 2]
rq2_cat2_col_id <- c()
rq2_cat1_col_id <- c()
for (col in 1:ncol(temp_1)) {
  idx1 <- which(!is.na(temp_1[,col]))
  # if 2 records from 2 different gas stations, then append compounds to cat2
  if (length(unique(transpose_mydata2_new[idx1]$gas_station)) > 1) {
    rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
  } else { # if 2 records from same gas stations, then append compounds to cat1
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
  # if 3 records from different gas stations, then append compounds to cat2
  if (length(unique(transpose_mydata2_new[idx2,]$gas_station)) > 1) {
    rq2_cat2_col_id <- c(rq2_cat2_col_id, col)
  } else {  # if 3 records from same gas stations, then append compounds to cat1
    rq2_cat1_col_id <- c(rq2_cat1_col_id, col)
  }
}

rq2_cat1 <- base::cbind(rq2_cat1, temp_2[, rq2_cat1_col_id])
rq2_cat2 <- base::cbind(rq2_cat2, temp_2[, rq2_cat2_col_id])

# if compounds have >=4 records -> it definitely appear in >=2 Gas stations -> append compounds to cat2
temp_3 <- transpose_mydata2_new[,3:ncol(transpose_mydata2_new)][,colSums(!is.na(transpose_mydata2_new[,3:ncol(transpose_mydata2_new)])) > 3]
rq2_cat2 <- base::cbind(rq2_cat2, temp_3)

# Check RQ2_category 1 (Compound found in only 1 gas station and not in any other) #####################################

# How many observations - 3046 obs
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% colnames(rq2_cat1)))

# What are the ASTM compounds?
ASTM_list[which(str_detect(ASTM_list$`RQ2: category 1: compounds with only 1 Gas record (rt1_thres = 0.2)`, "x")), "Compound"]

# "1,8-Dimethylnaphthalene" - Compound_8648.
# "1,2,3,5-Tetramethylbenzene" - Compound_3779.
# "2,3-Dimethylnaphthalene" - Compound_8125.
# "1-Ethylnaphthalene" - Compound_7736.
# "2-Ethylnaphthalene" - compound_7741.
# "Dodecane" - Compound_4510.
# "Undecane" - Compound_3336.
# "4-Isopropyltoluene" - Compound_2660.
# "3-Ethyltoluene" - Compound_1945.

# Check category 2 ( >=2 Gas Records) ##################################################################################
# How many observations - 11272 obs
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% colnames(rq2_cat2)))

# What are the ASTM compounds?
print(ASTM_list[which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 records and occur in >=2 Gas stations (1050 compounds) (rt1_thres = 0.2)`, "x")), "Compound"], n = 49)

# "Octane" 
# "Nonane"
# "Propylbenzene"
# "3-Ethyltoluene"
# "4-Ethyltoluene"
# "1,3,5-Trimethylbenzene"
# "2-Ethyltoluene"
# "Decane"
# "1,2,4-Trimethylbenzene"
# "Isobutylbenzene"
# "Sec-butylbenzene"
# "3-Isopropyltoluene"
# "4-Isopropyltoluene"
# "2-Isopropyltoluene"
# "Indane"
# "1,3-Diethylbenzene"
# "3- & 4-Propyltoluene"
# "1-Ethyl-3,5-Dimethylbenzene"
# "Undecane"
# "1-Ethyl-2,3-dimethylbenzene"
# "Isopentylbenzene"
# "Dodecane"
# "4,7-Dimethylindane"
# "Naphthalene"
# "Hexylbenzene"
# "Tridecane"
# "2-Methylnaphthalene"
# "1-Methylnaphthalene"
# "Tetradecane"
# "1,6-Dimethylnaphthalene"
# "2,3-Dimethylnaphthalene"
# "1,4-Dimethylnaphthalene"
# "2,6- & 1,3- & 1,7-Dimethylnaphthalene"
# "Isopropylbenzene"
# "1,2,3-Trimethylbenzene"
# " 2-Propyltoluene"
# "1-Ethyl-2,5-dimethylbenzene"
# "1-Ethyl-3,4-dimethylbenzene"
# "1-Ethyl-2,4-dimethylbenzene"
# "1-Ethyl-2,6-dimethylbenzene"
# "1,2,4,5-Tetramethylbenzene"
# "1-Methyl-2-tert-butylbenzene"
# "5-Methylindane"
# "2-Methylindane"
# "1,2,3,4-Tetramethylbenzene"
# "Pentylbenzene"
# "Dimethylindane (23.8366)"
# "Dimethylindane (24)"
# "Dimethylindane (24.1321)"

# Check Category 2 for presence of ASTM compounds BEFORE Wilcoxon test ###############################################################################################
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

# How many observations - 8292 peaks
View(shared_comp_normalized %>%
       filter(., collapsed_compound %in% colnames(rq2_cat2_stats)))

# How many ASTM compounds?
print(ASTM_list[which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION (506 compounds), BEFOREWilcoxon test (data imputed with LOD), alpha threshold < 0.1 (rt1_thres = 0.2)`, "x")), "Compound"], n = 49)


# Check category 2 for presence of ASTM compounds AFTER Wilcoxon test ###############################################################################################
# Alpha threshold = 0.05 - How many ASTM compounds? - 19/80
print(ASTM_list[which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION, AFTER Wilcoxon test (data imputed with LOD), alpha threshold < 0.05 (rt1_thres = 0.2)`, "x")), "Compound"], n = 49)

# Alpha threshold = 0.1 - How many ASTM compounds? - 23/80
print(ASTM_list[which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION, AFTER Wilcoxon test (data imputed with LOD), alpha threshold < 0.1 (rt1_thres = 0.2)`, "x")), "Compound"], n = 49)

# INtersect between PCA top loadings and compounds that passed Wilcoxon alpha = 0.1 ====================================
ASTM_intersect_rq2_pca_wilcox <- c()
for (ele in intersect_PCA_Wilcoxon_alpha0.1_rq2) {
  idx <- which(str_detect(ASTM_list$`RQ2: category 2: Compounds with >=2 Gas record FOR EACH GAS STATION, AFTER Wilcoxon test (data imputed with LOD), alpha threshold < 0.1 (rt1_thres = 0.2)`, 
                          ele))
  ASTM_intersect_rq2_pca_wilcox <- c(ASTM_intersect_rq2_pca_wilcox, 
                                     as.character(ASTM_list[idx, "Compound"]))
}

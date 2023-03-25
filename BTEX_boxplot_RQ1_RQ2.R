including <- function(df, include_list) {
  clean_data <- copy(df)
  temp <- list()
  i <- 1
  for (ele in c("^Benzene$",
                "^Toluene$",
                "^Ethylbenzene$",
                "Xylene")) {
    temp[[i]] <- clean_data %>%
      filter(grepl(ele, Compound))
    i <- i + 1
  }
  new_df <- dplyr::bind_rows(temp)
  return(new_df)
}

# DF1: Taking out BTEX compound with exact matching name form GC-MS
df_list <- list()
for (i in 1:length(df_list_step1.1)) {
  df_list[[i]] <- df_list_step1.1[[i]] %>%
    mutate(sample_name = ILR_file_list[[i]]) %>%
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

BTEX_only <- purrr::map(df_list, including, include_list = c("^Benzene$",
                                                             "^Toluene$",
                                                             "^Ethylbenzene$",
                                                             "Xylene")) 
  
BTEX_only_df <- dplyr::bind_rows(BTEX_only)


# DF2: Taking out BTEX compound with match RT1, RT2, Ion1, Ion2
df_2 <- purrr::map(df_list, filtering, filter_list = c("^Carbon disulfide$", 
                                                       "Cyclotrisiloxane..hexamethyl",
                                                       "Cyclotetrasiloxane..octamethyl",
                                                       "^Benzene$",
                                                       "^Toluene$",
                                                       "^Ethylbenzene$",
                                                       "Xylene")) 
df_2 <- dplyr::bind_rows(df_2)

BTEXASTM <- ASTM_withBTEX_F001A_aligned[c(1,2,4,5,7),]

temp <- list()
i <- 1
for (row in 1:nrow(BTEXASTM)) {
  rt1 <- BTEXASTM[row,]$RT1_aligned
  rt2 <- BTEXASTM[row,]$RT2_aligned
  ion1 <- BTEXASTM[row,]$Ion1
  ion2 <- BTEXASTM[row,]$Ion2
  
  
  # filter data by index
  idx_thres <- which(df_2$RT1 <= (rt1 + 0.1) & df_2$RT1 >= (rt1 - 0.1) & 
                       df_2$RT2 <= (rt2 + 0.1) & df_2$RT2 >= (rt2 - 0.1) & 
                       df_2$Ion1 <= (ion1 + 0.2) & df_2$Ion1 >= (ion1 - 0.2) & 
                       df_2$Ion2 <= (ion2 + 0.2) & df_2$Ion2 >= (ion2 - 0.2))
  
  temp[[i]] <- df_2[idx_thres,]
  i <- i + 1
}

df_2 <- dplyr::bind_rows(temp)

# Combine DF1 and DF2
full_BTEX <- base::rbind(BTEX_only_df, df_2)

# Boxplot RQ1: BTEX in Gas vs Diesel -----------------
BTEX_rq1 <- full_BTEX %>%
  filter(., fuel_type %in% c("Gas", "Diesel"))

colnames(BTEX_rq1)[colnames(BTEX_rq1) == 'fuel_type'] <- 'Fuel Type'

ggplot() +
  geom_boxplot(data = BTEX_rq1, aes(x=Compound, y = Area, colour = `Fuel Type`)) +
  theme_classic(base_size = 20)

# Boxplot RQ2: BTEX in Gas stations 1,3,8 vs Gas station 5,7,9 -----------------
BTEX_rq2 <- full_BTEX %>%
  filter(., fuel_type %in% c("Gas"))

colnames(BTEX_rq2)[colnames(BTEX_rq2) == 'gas_station'] <- 'Gas Station'

ggplot() +
  geom_boxplot(data = BTEX_rq2, aes(x=Compound, y = Area, colour = `Gas Station`)) +
  theme_classic(base_size = 20)

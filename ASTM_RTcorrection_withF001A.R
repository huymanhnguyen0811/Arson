# Cross check ASTM with F00187, F00187-2, F00187-3

setwd("C:/Users/huyng/Desktop/Huy Nguyen/PhD_EnSciMan_Ryerson_University/Arson project/Rproject/data")

file_list <- list.files(pattern = '*.xlsx')

# Pipe operator for isolating IL types
ILR_file_list <- file_list %>%
  .[!str_ends(., "check.xlsx")] %>%
  .[!str_detect(., "grouping_compounds")] %>%
  .[!str_detect(., "Station")] %>%
  .[!str_detect(., "Workflow")] %>%
  .[!str_detect(., "ASTM")] %>%
  .[!str_detect(., "data_table")] %>%
  .[str_detect(., "F00187")]

# Import IL samples to list
F001A <- purrr::map(ILR_file_list, read_xlsx,
                              sheet = "Results")

# df_step1.1 <- dplyr::bind_rows(df_list_step1.1)

# remove spaces in column names in df_list_step1.1
for (i in 1:length(F001A)) {
  colnames(F001A[[i]]) <- gsub(" ", "", colnames(F001A[[i]]))
}

step2 <- purrr::map(F001A, filtering, filter_list = c("^Carbon disulfide$", 
                                                                     "Cyclotrisiloxane..hexamethyl",
                                                                     "Cyclotetrasiloxane..octamethyl",
                                                                     "^Benzene$",
                                                                     "^Toluene$",
                                                                     "^Ethylbenzene$",
                                                                     "Xylene")) 

list_remaining_area <- limit_obser(step2, ILR_file_list, cap = 50000)[[1]]

step3 <- bind_rows(list_remaining_area) 

step4 <- grouping_comp_ver1(step3,
                            rt1thres = 0.2,
                            rt2thres = 0.125,
                            ion1thres = 0.05, # Ion 1 and 2 indicates molecular structure (2 most prevalent mass-to-charge)
                            ion2thres = 0.05)
# ASTM compound list
ASTM_list <- read_xlsx(paste0(getwd(), "/ASTM Compound List.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

# Filter ASTM light & mid range ----------------------
ASTM_light <- ASTM_list[, 1:5] %>%
  filter(., (RT1 < 20 & RT1 > 8 & RT2 < 4.4 & RT2 > 3.4)) 

ASTM_mid <- ASTM_list[, 1:5] %>%
  filter(., (RT1 < 30 & RT1 > 23 & RT2 < 5.4 & RT2 > 3.5)) 

ASTM_light_mid <- rbind(ASTM_light, ASTM_mid)

# Visual estimation of the adjustment
ggplotly(ggplot(data = step4,
                aes(x=RT1, y= RT2,
                    colour = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_light_mid, pch = 4,
                      aes(x=RT1, y= RT2, 
                          colour = "ASTM",
                          label = Compound)))

test <- step4 %>% 
  filter(., (RT1 < 20 & RT1 > 10 & RT2 < 4.5 & RT2 > 4))

# Adjusting ASTM Retention Time indices
ASTM_light_mid_F001A_aligned <- ASTM_light_mid %>%
  mutate(RT1_aligned = RT1-0.04) %>%
  mutate(RT2_aligned = RT2 -0.06)

# Visual inspection of the adjustment 
ggplotly(ggplot(data = step4,
                aes(x=RT1, y= RT2,
                    colour = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_light_mid_F001A_aligned, pch = 4,
                      aes(x=RT1_aligned, y= RT2_aligned, 
                          colour = "ASTM",
                          label = Compound)))

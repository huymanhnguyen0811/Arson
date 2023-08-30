# ASTM compound list
ASTM_list <- read_xlsx(paste0(getwd(), "/ASTM Compound List.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

# Filter ASTM light & mid range ----------------------
ASTM_light <- ASTM_list[, 1:5] %>%
  filter(., (RT1 < 20 & RT1 > 8 & RT2 < 4.4 & RT2 > 3.4)) 

ASTM_mid <- ASTM_list[, 1:5] %>%
  filter(., (RT1 < 30 & RT1 > 23 & RT2 < 5.4 & RT2 > 3.5)) 

ASTM_light_mid <- rbind(ASTM_light, ASTM_mid)

ASTM_F001A_aligned <- ASTM_list %>%
  mutate(RT1_aligned = RT1-0.04) %>%
  mutate(RT2_aligned = RT2 -0.06)

# Filter testing rt df light mid --------------
light_mid <- function(df) {
  testRT_light1 <- df %>%
    filter(., (RT1 < 20 & RT1 > 8 & RT2 < 3.6 & RT2 > 3.4))
  
  testRT_light2 <- df %>%
    filter(., (RT1 < 20 & RT1 > 8 & RT2 < 4.4 & RT2 > 4.1)) 
  
  testRT_mid1 <- df %>%
    filter(., (RT1 < 30 & RT1 > 23 & RT2 < 5.4 & RT2 > 5.1))
  
  testRT_mid2 <- df %>%
    filter(., (RT1 < 30 & RT1 > 23 & RT2 < 3.6 & RT2 > 3.5)) 
  
  testRT_mid3 <- df %>%
    filter(., (RT1 < 25 & RT1 > 23 & RT2 < 4.7 & RT2 > 4.1)) 
  
  testRT_light_mid <- rbind(testRT_light1, testRT_light2, testRT_mid1, testRT_mid2, testRT_mid3)
  return(testRT_light_mid)
}

test1_light_mid <- light_mid(shared_comp_normalized_rt10.1)
test2_light_mid <- light_mid(shared_comp_normalized_rt10.2)
test3_light_mid <- light_mid(shared_comp_normalized_rt10.3)
  
# Plots ----------------
  
ggplotly(ggplot(data = test1_light_mid,
                aes(x=RT1, y= RT2,
                    colour = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_F001A_aligned, pch = 4,
                      aes(x=RT1_aligned, y= RT2_aligned, 
                          colour = "ASTM",
                          label = Compound)))
           # scale_x_continuous(breaks = seq(min(test1_light_mid$RT1), max(test1_light_mid$RT1), 0.1)) +
           # scale_y_continuous(breaks = seq(min(test1_light_mid$RT2), max(test1_light_mid$RT2), 0.15)) +
           # theme(panel.grid = element_line(color = "#8ccde3",
           #                                 linewidth = 0.75,
           #                                 linetype = 2),
           #       axis.text.x = element_text(angle = 90)))

ggplotly(ggplot(data = test2_light_mid,
                aes(x=RT1, y= RT2,
                    colour = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_light_mid_F001A_aligned, pch = 4,
                      aes(x=RT1_aligned, y= RT2_aligned,
                          colour = "ASTM",
                          label = Compound)))
           # scale_x_continuous(breaks = seq(min(test1_light_mid$RT1), max(test1_light_mid$RT1), 0.2)) +
           # scale_y_continuous(breaks = seq(min(test1_light_mid$RT2), max(test1_light_mid$RT2), 0.15)) +
           # theme(panel.grid = element_line(color = "#8ccde3",
           #                                 linewidth = 0.75,
           #                                 linetype = 2),
           #       axis.text.x = element_text(angle = 90)))

ggplotly(ggplot(data = test3_light_mid,
                aes(x=RT1, y= RT2,
                    colour = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_light_mid_F001A_aligned, pch = 4,
                      aes(x=RT1_aligned, y= RT2_aligned,
                          colour = "ASTM",
                          label = Compound)))
           # scale_x_continuous(breaks = seq(min(test1_light_mid$RT1), max(test1_light_mid$RT1), 0.3)) +
           # scale_y_continuous(breaks = seq(min(test1_light_mid$RT2), max(test1_light_mid$RT2), 0.15)) +
           # theme(panel.grid = element_line(color = "#8ccde3",
           #                                 linewidth = 0.75,
           #                                 linetype = 2),
           #       axis.text.x = element_text(angle = 90)))
  

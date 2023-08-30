# ASTM compound list
ASTM_list <- read_xlsx(paste0(getwd(), "/ASTM Compound List.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

ASTM_list_withBTEX <- read_xlsx(paste0(getwd(), "/ASTM Compound List_withBTEX.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

ASTM_aligned_with_F001A <- ASTM_list %>%
  mutate(RT1_aligned = RT1 - 0.04) %>%
  mutate(RT2_aligned = RT2 - 0.06) %>%
  relocate(`RT1_aligned`, `RT2_aligned`, .after = RT2)

ASTM_withBTEX_aligned_with_F001A <- ASTM_list_withBTEX %>%
  mutate(RT1_aligned = RT1 - 0.04) %>%
  mutate(RT2_aligned = RT2 - 0.06) %>%
  relocate(`RT1_aligned`, `RT2_aligned`, .after = RT2)


# RQ1 ASTM confirm -------------------

ggplotly(ggplot(data = shared_comp_normalized_rt10.1 %>%
                  filter(., collapsed_compound %in% rq1_alpha0.1$collapsed_compound)) +
           geom_point(aes(x=RT1, y= RT2,
                          label = collapsed_compound),
                      pch = 1, size = 0.5) + 
           geom_point(data = ASTM_F001A_aligned, pch = 4, # ASTM_F001A_aligned
                      aes(x=RT1_aligned, y= RT2_aligned,
                          color = "ASTM",
                          label = Compound)) +
           theme(legend.position = "none") +
           theme_bw(base_size = 20) +
           xlim(5,36))



# RQ2 ASTM confirm ----------------
ggplotly(ggplot(data = shared_comp_normalized_rt10.1 %>%
                  filter(., collapsed_compound %in% wilcox_result_rt10.1[[2]]$collapsed_compound),
                aes(x=RT1, y= RT2,
                    label = collapsed_compound)) +
           geom_point(pch = 1, size = 0.5) + 
           geom_point(data = ASTM_withBTEX_F001A_aligned, pch = 4, # ASTM_F001A_aligned
                      aes(x=RT1_aligned, y= RT2_aligned,
                          color = "ASTM Compounds",
                          label = Compound)))


# ggplot(data = ASTM_F001A_aligned, pch = 4,
#        aes(x=RT1_aligned, y= RT2_aligned)) +
#   # ggrepel::geom_label_repel(aes(label = Compound),
#   #                           max.overlaps = 100,
#   #                           force = 10) +
#   geom_point(pch = 4, colour = "red") + 
#   geom_point(data = shared_comp_normalized_rt10.2 %>%
#                filter(., collapsed_compound %in% wilcox_result_rt10.2[[2]]$collapsed_compound),
#              aes(x=RT1, y= RT2), size = 1, colour = "black")

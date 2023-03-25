# ASTM compound list
ASTM_list <- read_xlsx(paste0(getwd(), "/ASTM Compound List.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

ASTM_list_withBTEX <- read_xlsx(paste0(getwd(), "/ASTM Compound List_withBTEX.xlsx"), sheet = "Sheet1") %>%
  arrange(RT1, RT2)

ASTM_F001A_aligned <- ASTM_list %>%
  mutate(RT1_aligned = RT1 - 0.04) %>%
  mutate(RT2_aligned = RT2 - 0.06) %>%
  relocate(`RT1_aligned`, `RT2_aligned`, .after = RT2)

ASTM_withBTEX_F001A_aligned <- ASTM_list_withBTEX %>%
  mutate(RT1_aligned = RT1 - 0.04) %>%
  mutate(RT2_aligned = RT2 - 0.06) %>%
  relocate(`RT1_aligned`, `RT2_aligned`, .after = RT2)


# RQ1 ASTM confirm -------------------

ggplotly(ggplot(data = shared_comp_normalized_rt10.1 %>%
                  filter(., collapsed_compound %in% rq1_alpha0.1$collapsed_compound)) +
           geom_point(aes(x=RT1, y= RT2,
                          label = collapsed_compound),
                      pch = 1, size = 0.5) + 
           geom_point(data = ASTM_withBTEX_F001A_aligned, pch = 4, # ASTM_F001A_aligned
                      aes(x=RT1_aligned, y= RT2_aligned,
                          color = "ASTM Compounds",
                          label = Compound)))

beyondASTMrq1 <- rq1_alpha0.1 %>% filter(., collapsed_compound %notin% c("Compound_481.", "Compound_2269.", "Compound_4588.",
                                                                         "Compound_1557.", "Compound_1893.", "Compound_1981.",
                                                                         "Compound_2005.", "Compound_2065.", "Compound_2166.", 
                                                                         "Compound_2377.", "Compound_2758.", "Compound_3011.",
                                                                         "Compound_3045.", "Compound_3454.", "Compound_4259.",
                                                                         "Compound_4664.", "Compound_4823.", "Compound_4781.",
                                                                         "Compound_2969.", "Compound_4862.", "Compound_7863.", 
                                                                         "Compound_8048.", "Compound_8261.", "Compound_8327.",
                                                                         "Compound_8581."))
beyondASTMrq1_df <- shared_comp_normalized_rt10.1 %>% 
  filter(., collapsed_compound %in% unique(beyondASTMrq1$collapsed_compound)) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = median(RT1), 
            rt2 = median(RT2),
            ion1 = median(Ion1),
            ion2 = median(Ion2)) %>%
  arrange(rt1, rt2)

colnames(beyondASTMrq1_df) <- c("Collapsed_compound", "Retention time 1", "Retention time 2", 
                                "Molecular Ion1", "Molecular Ion2")

writexl::write_xlsx(beyondASTMrq1_df, path = paste0(getwd(), "/beyondASTMrq1_df.xlsx"))

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

beyondASTMrq2 <- wilcox_result_rt10.1[[2]] %>% filter(., collapsed_compound %notin% c("Compound_2269.", "Compound_1557.", "Compound_1893.",
                                                                                      "Compound_1981.", "Compound_2005.", "Compound_2065.",
                                                                                      "Compound_2166.", "Compound_2377.", "Compound_2508.",
                                                                                      "Compound_2680.", "Compound_2754.", "Compound_2758.",
                                                                                      "Compound_3011.", "Compound_3045.", "Compound_3129.",
                                                                                      "Compound_3233.", "Compound_3334.", "Compound_3373.",
                                                                                      "Compound_3504.", "Compound_3454.", "Compound_3691.",
                                                                                      "Compound_3803.", "Compound_3863.", "Compound_4259.", 
                                                                                      "Compound_4245.", "Compound_4722.", "Compound_4823.", 
                                                                                      "Compound_2969.", "Compound_4862.", "Compound_8261.", 
                                                                                      "Compound_8327."))

beyondASTMrq2_df <- shared_comp_normalized_rt10.1 %>% 
  filter(., collapsed_compound %in% unique(beyondASTMrq2$collapsed_compound)) %>%
  group_by(collapsed_compound) %>%
  summarise(rt1 = median(RT1), 
            rt2 = median(RT2),
            ion1 = median(Ion1),
            ion2 = median(Ion2)) %>%
  arrange(rt1, rt2)

colnames(beyondASTMrq2_df) <- c("Collapsed_compound", "Retention time 1", "Retention time 2", 
                                "Molecular Ion1", "Molecular Ion2")

writexl::write_xlsx(beyondASTMrq2_df, path = paste0(getwd(), "/beyondASTMrq2_df.xlsx"))


# ggplot(data = ASTM_F001A_aligned, pch = 4,
#        aes(x=RT1_aligned, y= RT2_aligned)) +
#   # ggrepel::geom_label_repel(aes(label = Compound),
#   #                           max.overlaps = 100,
#   #                           force = 10) +
#   geom_point(pch = 4, colour = "red") + 
#   geom_point(data = shared_comp_normalized_rt10.2 %>%
#                filter(., collapsed_compound %in% wilcox_result_rt10.2[[2]]$collapsed_compound),
#              aes(x=RT1, y= RT2), size = 1, colour = "black")

library(lsa)
library(pheatmap)

# data matrix with Analytes are row wise, while samples are column wise. MUST NOT INCLUDE NA values

cos_theta_df <- base::cbind(rq2_cat2_stats, transpose_mydata2_new[, 1:2]) %>%
  relocate(sample_name, gas_station, .before = everything()) %>%
  select(., -gas_station) %>% 
  column_to_rownames(., var = "sample_name")

transpose_cos_theta_df <- data.table::transpose(cos_theta_df)
rownames(transpose_cos_theta_df) <- colnames(cos_theta_df)
colnames(transpose_cos_theta_df) <- rownames(cos_theta_df)

temp1 <- as.matrix(transpose_cos_theta_df)
temp2 <- as.matrix(lsa::cosine(temp1))


pheatmap(temp2, display_numbers = T, fontsize_number = 12, fontsize = 12)
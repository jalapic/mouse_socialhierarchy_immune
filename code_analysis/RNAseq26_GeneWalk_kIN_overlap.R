
`%notin%` <- Negate(`%in%`)

my_cutoff = 0.15


kIN_df_all <- readRDS(glue("results_RNAseqRDS/kIN_dataframe_Liver.RDS")) %>%
  mutate(tissue = "Liver") %>% 
  rbind(readRDS(glue("results_RNAseqRDS/kIN_dataframe_Spleen.RDS")) %>%
          mutate(tissue = "Spleen")) %>% 
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(desc(kWithin)))


kIN_df_all %>% 
  count(tissue, module, sort = T) 

genewalk_result <- readRDS("results_RNAseqRDS/genewalk_result.RDS")


table(genewalk_result$gene_group,genewalk_result$result_type)


genewalk_result %>% 
  select(tissue, module) %>% 
  unique() %>% 
  arrange(tissue,module) %>% 
  filter(module %notin% c("UP","DOWN"))-> temp_df

temp_df
i = 2

key_reg_list <- vector('list', length = nrow(temp_df))

for(i in 1:nrow(temp_df)){
  my_tissue = temp_df$tissue[i]
  my_module = temp_df$module[i]
  
  genewalk_reg <- genewalk_result %>% 
    filter(result_type == "Regulator") %>% 
    filter(tissue == my_tissue) %>% 
    filter(module == my_module) %>% 
    .$symbol
  
  kIN_df <- kIN_df_all %>% 
    filter(tissue == my_tissue) %>% 
    filter(module == my_module) %>% 
    arrange(desc(kWithin)) %>% 
    mutate(kINrank = rank(desc(kWithin)))
  
  
  genewalk_reg 
  
  kIN_df %>% 
    unique() %>% 
  nrow()
  
  round(nrow(kIN_df)*my_cutoff) -> cutoff
  
  
  genewalk_reg -> x
  
  kIN_df %>% 
    head(cutoff) %>% 
    .$symbol -> y
  
  x[x %in% y] -> overlap
  overlap
  
  kIN_df %>% 
    head(cutoff) %>%
    filter(symbol %in% overlap) %>% 
    filter(module == my_module) %>% 
    # select(module,symbol) %>% 
    unique() %>% 
    mutate(tissue = my_tissue) -> key_reg_list[[i]]
  
}


key_reg_list %>% 
  rbindlist %>% 
  write.csv(glue("results_tables/kIN_genewalk_regulator_ovelap_cutoff{my_cutoff}.csv"),row.names = F)



# moonlighters

genewalk_result %>% 
  filter(result_type == "Moonlighter") %>% 
  write.csv("results_tables/genewalk_moonlighters.csv",row.names = F)


genewalk_result %>% 
  filter(module %notin% c("UP","DOWN")) %>% 
  filter(result_type == "Moonlighter") %>% 
  count(gene_group)

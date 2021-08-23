kIN_df_all <- readRDS(glue("results_RNAseqRDS/kIN_dataframe_Liver.RDS")) %>%
  mutate(tissue = "Liver") %>% 
  rbind(readRDS(glue("results_RNAseqRDS/kIN_dataframe_Spleen.RDS")) %>%
          mutate(tissue = "Spleen")) %>% 
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(desc(kWithin)))


kIN_df_all %>% 
  filter(module %in% c("pink","cyan","blue","brown","green")) %>% 
  split(.$module) %>% 
  lapply(., function (x) x$symbol)-> symbol_list


limma_list_liver <- readRDS("results_RNAseqRDS/limma_Liver.RDS")
limma_list_spleen <- readRDS("results_RNAseqRDS/limma_Spleen.RDS")


y_liver <- readRDS("results_RNAseqRDS/limma_vdl_Liver")
y_spleen <- readRDS("results_RNAseqRDS/limma_vdl_Spleen")

i = 1

wgcna_deg_status <- vector('list',length(symbol_list))

for (i in 1:length(symbol_list)){
  
  names(symbol_list)[[i]] -> names(wgcna_deg_status)[[i]] -> ma_mod
  symbol_list[[i]] 
  
  if(ma_mod %in% c("pink","cyan",'blue')){
    limma_list_liver$status -> deg_df
  } else {limma_list_spleen$status -> deg_df}
  
  deg_df %>% 
    filter(symbol %in% symbol_list[[i]]) %>% 
    select(symbol,logFC, P.Value) %>% 
    arrange(desc(abs(logFC))) %>% 
    filter(P.Value <0.05) %>% 
    unique() %>% 
    head(5) -> wgcna_deg_status[[i]]
  
}

wgcna_deg_status %>% 
  map2_df(., names(.), ~mutate(.x, module = .y)) %>% 
  mutate(module = factor(module, levels = c("pink","cyan","blue","brown","green"))) %>% 
  arrange(module) %>% 
  mutate(DEG_by = "Status") -> deg_status



wgcna_deg_cort <- vector('list',length(symbol_list))

for (i in 1:length(symbol_list)){
  
  names(symbol_list)[[i]] -> names(wgcna_deg_cort)[[i]] -> ma_mod
  symbol_list[[i]] 
  
  if(ma_mod %in% c("pink","cyan",'blue')){
    limma_list_liver$cort -> deg_df
  } else {limma_list_spleen$cort -> deg_df}
  
  deg_df %>% 
    filter(symbol %in% symbol_list[[i]]) %>% 
    select(symbol,logFC, P.Value) %>% 
    arrange(desc(abs(logFC))) %>% 
    filter(P.Value <0.05) %>% 
    unique() %>% 
    head(5) -> wgcna_deg_cort[[i]]
  
}

wgcna_deg_cort %>% 
  map2_df(., names(.), ~mutate(.x, module = .y)) %>% 
  mutate(module = factor(module, levels = c("pink","cyan","blue","brown","green"))) %>% 
  arrange(module) %>% 
  mutate(DEG_by = "CORT") -> deg_cort



rbind(deg_status,deg_cort) -> deg_df

write.csv(deg_df, "results_tables/wgcna_top_DEG.csv", row.names = F)

deg_df$symbol -> my_genes



deg_status$symbol -> my_genes
grcm38 %>% 
  filter(symbol %in% my_genes) %>% 
  .$ensgene %>% 
  unique() -> gois

y_liver$E[gois,] %>% t %>% 
  as.data.frame() %>%   rownames_to_column("subjectID") %>% 
  gather(ensgene,Expression,2:(ncol(.))) %>% 
  left_join(alldata %>% 
              select(subjectID,domgroup)) %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) %>% 
  left_join(deg_status) -> deg_plot_df

deg_plot_df %>% filter(is.na(logFC))

deg_plot_df %>% 
  filter(as.numeric(Expression) > 20)

hist(deg_plot_df$Expression)

deg_plot_df %>% 
  ggplot(aes(domgroup,reorder(symbol,abs(logFC)), fill = Expression))+
  geom_tile()+
  scale_fill_distiller(palette = "PiYG") +
  facet_grid(module ~., scales = "free_y")+
  theme_minimal()
  







library(pheatmap)

y_liver$E[gois,] %>% t  
  
  
  
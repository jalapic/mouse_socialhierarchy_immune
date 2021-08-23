genewalk_names <- list.files(path ="results_GeneWalk",pattern = 'immune_*')
genewalk_shortnames <- gsub("immune_","",genewalk_names)


genewalkresults_list <- lapply(glue('results_GeneWalk/{genewalk_names}/genewalk_results.csv'), 
                               function (x) read_csv(x) %>% 
                                 select(hgnc_symbol,entrez_mouse) %>% 
                                 rename(entrez = entrez_mouse))


moonlighter_list <- lapply(glue('results_GeneWalk/{genewalk_names}/figures/genewalk_moonlighters.csv'), read_csv)

regulator_list <- lapply(glue('results_GeneWalk/{genewalk_names}/figures/genewalk_regulators.csv'), read_csv)

names(genewalkresults_list) <- genewalk_shortnames
names(moonlighter_list) <- genewalk_shortnames
names(regulator_list) <- genewalk_shortnames


i = 1 
for (i in 1:length(genewalk_names)){
  moonlighter_list[[i]] %>% 
    rename(hgnc_symbol = gw_moonlighter) %>% 
    left_join(genewalkresults_list[[i]]) %>% 
    left_join(grcm38 %>% select(entrez, symbol)) %>% 
    unique() -> moonlighter_list[[i]]
}
  
for (i in 1:length(genewalk_names)){
  regulator_list[[i]] %>% 
    rename(hgnc_symbol = gw_regulators) %>% 
    left_join(genewalkresults_list[[i]]) %>% 
    left_join(grcm38 %>% select(entrez, symbol)) %>% 
    unique() -> regulator_list[[i]]
}


moonlighter_list %>% 
  map2_df(.,names(.),~mutate(.x, gene_group =.y )) %>% 
  select(symbol, gene_group) %>% 
  mutate(result_type = "Moonlighter") -> x 

regulator_list %>% 
  map2_df(.,names(.),~mutate(.x, gene_group =.y )) %>% 
  select(symbol, gene_group) %>% 
  mutate(result_type = "Regulator") -> y

rbind(x, y) %>% 
  mutate(module = gsub(".*_", "", gene_group)) %>% 
  mutate(tissue = gsub("\\_.*", "", gene_group)) -> genewalk_result

saveRDS(genewalk_result, "results_RNAseqRDS/genewalk_result.RDS")


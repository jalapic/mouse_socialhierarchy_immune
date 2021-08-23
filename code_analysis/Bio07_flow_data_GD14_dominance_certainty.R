#take a look at dominance certainty ~ rank first 



# now flow data 
flow_GD14 <- flow_all %>% 
  filter(housing=="group") %>% 
  filter(cohort!="Pair-comparison") %>% 
  left_join(alldata %>% select(subjectID, ds, glicko_rank, ind_cert)) %>% 
  mutate(glicko_rank=as.integer(glicko_rank))


#stat

list_flow_GD14 <- flow_GD14 %>% 
  split(.$cell_type) %>% 
  map(~ mutate(., perc = scale(perc)))

result_list_domcert_GD14 <- list()

for (i in 1:length(list_flow_GD14)){
  
  result_list_domcert_GD14[[i]] <- brm(
    perc ~ ind_cert + (1|cohort),
    file = glue ('results_statRDS/GD14_domcert_{names(list_flow_GD14)[[i]]}.Rds'),
    data = list_flow_GD14[[i]])
}

names(result_list_domcert_GD14) <- names(list_flow_GD14)



# main manuscript figure ============================================

result_list_domcert_GD14[[1]] -> m
get_variables(m)

lapply(result_list_domcert_GD14, 
       function(x) x %>% 
         spread_draws(b_ind_cert) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_gd14_domcert


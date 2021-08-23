%>% %>% flow_PD09<-flow_all %>% 
  filter(housing=="pair") %>% 
  filter(!is.na(pair_status)) %>% 
  mutate(dyadID=paste(batch,mouseID,sep=".")) %>% 
  left_join(alldata %>% select(subjectID, ds, glicko_rank)) %>% 
  mutate(glicko_rank=as.integer(glicko_rank))


list_flow_PD09<-flow_PD09 %>% 
  split(.$cell_type)%>% 
  map(~ mutate(., perc = scale(perc)))

result_list_predict <- list()
for (i in 1:length(list_flow_PD09)){
  
  result_list_predict[[i]] <- brm(
    ds~perc + (1|cohort),
    family = gaussian,
    file = glue ('results_statRDS/PD09_predict_{names(list_flow_PD09)[[i]]}.Rds'),
    data = list_flow_PD09[[i]])
}
names(result_list_predict) <- names(list_flow_PD09)


# main manuscript figure

result_list_predict[[1]] -> m
get_variables(m)

lapply(result_list_predict, 
       function(x) x %>% 
         spread_draws(b_perc) %>%
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_pd09_predict




result_list_predict_glick[[1]] -> m
get_variables(m)

lapply(result_list_predict_glick, 
       function(x) x %>% 
         spread_draws(b_perc) %>%
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_pd09_predict_glick


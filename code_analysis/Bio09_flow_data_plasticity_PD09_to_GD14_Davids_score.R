colnames(alldata)

flow_all %>% 
  filter(cohort!="Pair-comparison") %>% 
  select(subjectID,housing,cell_type,perc) %>% 
  left_join(alldata %>% select(subjectID,domgroup, ds)) %>% 
  mutate(domgroup = factor(domgroup, levels = c("Alpha","Subdominant", "Subordinate"),ordered = T)) %>% 
  mutate(housing = factor(housing, levels = c("pair","group"))) %>% 
  split(.$cell_type) %>% 
  map(~ mutate(., perc_raw = perc)) %>% 
  map(~ mutate(., perc = scale(perc)))  -> list_flow_plasticity_2

levels(factor(list))
# David'score
result_list_plasticity_ds <- list()

for (i in 1:length(list_flow_plasticity_2)){
  
  result_list_plasticity_ds[[i]] <- brm(
    perc ~ ds+housing+(housing||subjectID),
    family = gaussian,
    file = glue ('results_statRDS/plasticity_ds_{names(list_flow_plasticity_2)[[i]]}.Rds'),
    data = list_flow_plasticity_2[[i]])
}

names(result_list_plasticity_ds) <- names(list_flow_plasticity_2)


# domgroup: social status group 
result_list_plasticity_domgroup <- list()

for (i in 1:length(list_flow_plasticity_2)){
  
  result_list_plasticity_domgroup[[i]] <- brm(
    perc ~ mo(domgroup)+housing+(housing||subjectID),
    family = gaussian,
    file = glue ('results_statRDS/plasticity_domgroup_{names(list_flow_plasticity_2)[[i]]}.Rds'),
    data = list_flow_plasticity_2[[i]])
}

names(result_list_plasticity_domgroup) <- names(list_flow_plasticity_2)



flow_GD14 <- flow_all %>% 
  filter(housing=="group") %>% 
  filter(cohort!="Pair-comparison") %>% 
  left_join(alldata %>% select(subjectID, ds, glicko_rank, Rating)) %>% 
  mutate(glicko_rank=as.integer(glicko_rank))



#stat

list_flow_GD14 <- flow_GD14 %>% 
  split(.$cell_type)%>% 
  map(~ mutate(., perc = scale(perc)))

result_list_ds_GD14 <- list()

for (i in 1:length(list_flow_GD14)){
  
  result_list_ds_GD14[[i]] <- brm(
    perc ~ ds + (1|cohort),
    file = glue ('results_statRDS/GD14_ds_{names(list_flow_GD14)[[i]]}.Rds'),
    data = list_flow_GD14[[i]])
}

names(result_list_ds_GD14) <- names(list_flow_GD14)

result_list_ds_GD14[[1]] -> m
get_variables(m)

lapply(result_list_ds_GD14, 
       function(x) x %>% 
         spread_draws(b_ds) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_gd14_ds




# stat - monotonic effect - rank 


list_flow_GD14 <- flow_GD14 %>% 
  split(.$cell_type)%>% 
  map(~ mutate(., perc = scale(perc)))

result_list_moglick_GD14 <- list()

for (i in 1:length(list_flow_GD14)){
  
  result_list_moglick_GD14[[i]] <- brm(
    perc ~ mo(glicko_rank) + (1|cohort),
    file = glue ('results_statRDS/GD14_mo_glickorank_{names(list_flow_GD14)[[i]]}.Rds'),
    data = list_flow_GD14[[i]])
}

names(result_list_moglick_GD14) <- names(list_flow_GD14)

result_list_moglick_GD14[[3]] -> m
get_variables(m)

lapply(result_list_moglick_GD14, 
       function(x) x %>% 
         spread_draws(bsp_moglicko_rank) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_gd14_moglick



conditional_effects(m)



#stat - glicko rating instead of David's score 

list_flow_GD14 <- flow_GD14 %>% 
  split(.$cell_type)%>% 
  map(~ mutate(., perc = scale(perc))) %>% 
  map(~ mutate(., scaled_Rating = scale(Rating)))

result_list_glickorating_GD14 <- list()

for (i in 1:length(list_flow_GD14)){
  
  result_list_glickorating_GD14[[i]] <- brm(
    perc ~ scaled_Rating + (1|cohort),
    file = glue ('results_statRDS/GD14_glickorating_{names(list_flow_GD14)[[i]]}.Rds'),
    data = list_flow_GD14[[i]])
}

names(result_list_glickorating_GD14) <- names(list_flow_GD14)

result_list_glickorating_GD14[[1]] -> m
get_variables(m)

lapply(result_list_glickorating_GD14, 
       function(x) x %>% 
         spread_draws(b_scaled_Rating) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_gd14_glickorating



#sUBORDINATE ONLY: MONOCYTE
alldata %>% 
  filter(domgroup == "Subordinate") -> sub



png(filename = glue("results_figures/monogd14_lossprop_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(loss_prop,Mono_GD14, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  scale_fill_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Loss proportion", 
       y = "Monocyte proportion (GD14)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())


png(filename = glue("results_figures/monogd14_winprop_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(win_prop,Mono_GD14, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  scale_fill_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Win proportion", 
       y = "Monocyte proportion (GD14)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())


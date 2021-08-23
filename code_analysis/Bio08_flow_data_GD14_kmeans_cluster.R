
colnames(alldata)



flow_GD14 <- flow_all %>% 
  filter(housing=="group") %>% 
  filter(cohort!="Pair-comparison") %>% 
  left_join(alldata %>% select(subjectID, km3_cluster, domgroup,pair_status)) 



#stat

list_flow_GD14 <- flow_GD14 %>% 
  split(.$cell_type)%>% 
  map(~ mutate(., perc = scale(perc)))

result_list_cluster_GD14 <- list()

for (i in 1:length(list_flow_GD14)){
  
  result_list_cluster_GD14[[i]] <- brm(
    perc ~ km3_cluster,
    file = glue ('results_statRDS/GD14_cluster_{names(list_flow_GD14)[[i]]}.Rds'),
    data = list_flow_GD14[[i]])
}

names(result_list_cluster_GD14) <- names(list_flow_GD14)

result_list_cluster_GD14[[1]] -> m
get_variables(m)

lapply(result_list_cluster_GD14, 
       function(x) x %>% 
         spread_draws(b_km3_clusterCluster2) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_gd14_cluster



# plot =================================================================
mycellname <- data.frame (
  cell_type =  names(list_flow_GD14),
  pretty_cell_type = c("B cells", "CD4/CD8", "Cytotoxic T", "DCs",
                       "Helper T", "Lymphoid DCs", "Macrophages","MLR",
                       "Monocytes", "Myeloid DCs","Neutrophils","NK cells",
                       "NLR","T cells")) %>% 
  mutate(pretty_cell_type = factor(pretty_cell_type,
                                   levels = c("Macrophages","Monocytes","Neutrophils",
                                                  "NK cells","B cells","T cells",
                                                  "Helper T","Cytotoxic T",
                                                  "Lymphoid DCs","Myeloid DCs", "DCs",
                                                  "CD4/CD8","MLR","NLR")))



flow_GD14 %>% 
  left_join(mycellname) %>% 
  mutate(myx = domgroup) %>% 
  filter(pretty_cell_type!="DCs") %>% 
    mutate(myx = km3_cluster) %>% 
  ggplot(aes(myx,perc,color=myx,fill=myx))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 jitter.height = 0.01, jitter.width = 0.1, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  facet_wrap( ~ pretty_cell_type, scales = "free_y", nrow = 3)+
  labs(x = "Cluster",
       y = "% of cells",
       color = "Cluster",
       fill = "Cluster")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none") -> cluster_GD14


flow_GD14 %>% 
  left_join(mycellname) %>% 
  mutate(myx = domgroup) %>% 
  filter(pretty_cell_type!="DCs") %>% 
  ggplot(aes(myx,perc,color=myx,fill=myx))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 jitter.height = 0.01, jitter.width = 0.1, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  facet_wrap( ~ pretty_cell_type, scales = "free_y", nrow = 3)+
  labs(x = "Cluster",
       y = "% of cells",
       color = "Social status group",
       fill = "Social status group")+
  scale_color_viridis(discrete = T, option = "D")+
  scale_fill_viridis(discrete = T, option = "D")+
  theme_bw()+
  theme(legend.position = "none") -> domgroup_GD14


ggsave(cluster_GD14, filename = "results_figures/cluster_GD14.png",width = 12,height = 7.5)
ggsave(domgroup_GD14, filename = "results_figures/domgroup_GD14.png",width = 12,height = 7.5)

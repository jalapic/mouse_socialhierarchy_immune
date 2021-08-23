
colnames(flow_all)

flow_PD09<-flow_all %>% 
  filter(housing=="pair") %>% 
  filter(!is.na(pair_status)) %>% 
  mutate(dyadID=paste(batch,mouseID,sep=".")) %>% 
  left_join(alldata %>% select(subjectID, ds, glicko_rank)) %>% 
  mutate(glicko_rank=as.integer(glicko_rank))


#Robust estmiation with brms ===================================================
table(flow_PD09$cell_type)

list_flow_PD09<-flow_PD09 %>% 
  split(.$cell_type) %>% 
  map(~ mutate(., perc = scale(perc)))


result_list <- list()

for (i in 1:length(list_flow_PD09)){
  
  result_list[[i]] <- brm(
    bf(perc~pair_status,
       sigma ~ pair_status),
    # control = list(adapt_delta=0.95),
    family = student,
    file = glue ('results_statRDS/PD09_pairstatus_{names(list_flow_PD09)[[i]]}.Rds'),
    data = list_flow_PD09[[i]])
}
names(result_list) <- names(list_flow_PD09)


# main manuscript figure

result_list[[1]] -> m
get_variables(m)

lapply(result_list, 
       function(x) x %>% 
         spread_draws(b_pair_statussub) %>%
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_pd09_pairstatus

   


# Just do mixed effect model with dom-sub category===================================================
table(flow_PD09$cell_type)

list_flow_PD09<-flow_PD09 %>% split(.$cell_type)
result_list <- list()


for (i in 1:length(list_flow_PD09)){

  result_list[[i]] <- brm(
  bf(perc~pair_status + (1|dyadID)),
  # control = list(adapt_delta=0.95),
  family = gaussian,
  file = glue ('results_statRDS/PD09_pairstatus_{names(list_flow_PD09)[[i]]}.Rds'),
  data = list_flow_PD09[[i]])
  }
names(result_list) <- names(list_flow_PD09)

# now extract estimates
pairstatus_PD09_result = data.frame(cell_type = names(list_flow_PD09),Estimate = NA,lower = NA, upper = NA)

for (i in 1:length(list_flow_PD09)){
  result_list[[i]] %>% summary() -> x
  x$fixed %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname == "pair_statussub") %>%
    select(Estimate,`l-95% CI`,`u-95% CI`) -> pairstatus_PD09_result[i,2:4]
}

pairstatus_PD09_result %>%
  mutate(sig = ifelse(lower*upper >=0, "Sig","")) %>% 
  mutate(stat = glue('{round(Estimate,2)} [{round(lower,2)}, {round(upper,2)}]')) -> pairstatus_PD09_result


# PLOT ============================================================================
list_dat<-list()
list_plot<-list()

for(i in 1:14){
  list_dat[[i]]<-conditional_effects(result_list[[i]])$pair_status %>%
    mutate(xvar=pair_status,yvar=estimate__,lower=lower__,upper=upper__) %>% 
    select(xvar,yvar,lower,upper) %>% 
    as.data.frame()
  
  list_plot[[i]]<-ggplot(list_flow_PD09[[i]],aes(pair_status,perc))+
    geom_line(aes(group=dyadID),color="grey80",size=0.8)+
    geom_point(aes(color=pair_status,fill=pair_status),shape=21,size=4.8,alpha=0.4)+
    labs(x="",
         y="")+
    scale_x_discrete(labels=c("Dom","Sub"))+
    scale_color_viridis(discrete = T)+
    scale_fill_viridis(discrete = T)+
    annotate("pointrange",
             x = 2.2 ,
             y = list_dat[[i]] %>% filter(xvar=="sub") %>% .$yvar, 
             ymin = list_dat[[i]] %>% filter(xvar=="sub") %>% .$lower, 
             ymax = list_dat[[i]] %>% filter(xvar=="sub") %>% .$upper,
             size = 1.5,
             color="#FDE725FF",alpha=0.9)+
    annotate("pointrange",
             x = 0.8, 
             y = list_dat[[i]] %>% filter(xvar=="dom") %>% .$yvar, 
             ymin = list_dat[[i]] %>% filter(xvar=="dom") %>% .$lower, 
             ymax = list_dat[[i]] %>% filter(xvar=="dom") %>% .$upper,
             size = 1.5,
             color="#440154FF",alpha=0.9)+
    ggtitle(glue('{names(result_list)[[i]]}'))+
    theme_base(base_size = 15)+
    theme(legend.position = "none",
          plot.title = element_text(size = 13),
          plot.background = element_blank()
          )
  
  
}



pairplot <- grid.arrange(list_plot[[1]] + ggtitle(glue('B cells {pairstatus_PD09_result$stat[1]}')),
             list_plot[[2]] + ggtitle(glue('CD4/CD8 ratio {pairstatus_PD09_result$stat[2]}')),
             list_plot[[3]] + ggtitle(glue('cytotoxic T {pairstatus_PD09_result$stat[3]}')),
             list_plot[[4]] + ggtitle(glue('Dendritic cells {pairstatus_PD09_result$stat[4]}')),
             list_plot[[5]] + ggtitle(glue('Helper T {pairstatus_PD09_result$stat[5]}')),
             list_plot[[6]] + ggtitle(glue('Lymphoid DCs {pairstatus_PD09_result$stat[6]}')),
             list_plot[[7]] + ggtitle(glue('Marcophages {pairstatus_PD09_result$stat[7]}')),
             list_plot[[8]] + ggtitle(glue('MLR {pairstatus_PD09_result$stat[8]}')),
             list_plot[[9]] + ggtitle(glue('Monocytes {pairstatus_PD09_result$stat[9]}')),
             list_plot[[10]] + ggtitle(glue('Myeloid DCs {pairstatus_PD09_result$stat[10]}')),
             list_plot[[11]] + ggtitle(glue('Neutrophils {pairstatus_PD09_result$stat[11]}')),
             list_plot[[12]] + ggtitle(glue('NK cells {pairstatus_PD09_result$stat[12]}')),
             list_plot[[13]] + ggtitle(glue('NLR {pairstatus_PD09_result$stat[13]}')),
             list_plot[[14]] + ggtitle(glue('T cells {pairstatus_PD09_result$stat[14]}')),
             nrow = 2)

ggsave(pairplot, filename = "results_figures/pairstatus_PD09.png",width = 25,height = 11) # nrow = 2


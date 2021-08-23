# make sure you run Bio09_flow_data_plasticity_PD09_to_GD14.R first. 

flow_all %>% 
  filter(cohort!="Pair-comparison") %>% 
  select(subjectID,housing,cell_type,perc) %>%
  left_join(alldata %>% select(subjectID,domgroup,ds)) %>% 
  mutate(housing = factor(housing, levels = c("pair","group"))) %>% 
  split(.$cell_type)-> list_dat

names(list_dat) <- c("B cells", "CD4/CD8", "Cytotoxic T", "DCs",
  "Helper T", "Lymphoid DCs", "Macrophages","MLR",
  "Monocytes", "Myeloid DCs","Neutrophils","NK cells",
  "NLR","T cells")


get_coordinate <- function(df){
  conditional_effects(df)$housing %>% 
    select(housing,lower__,upper__) -> x
  x$estimate__ <-conditional_effects(df)$housing$`estimate__`
  colnames(x) <- c("housing", "lower", "upper", "estimate")
  return(x)
}

list() -> cond_dat
list() -> plasticity_plot_list_domgroup
list() -> plasticity_plot_list

# i = 14
# result_list_plasticity_domgroup[[i]] -> df


for (i in 1:14){
cond_dat[[i]] <- get_coordinate(result_list_plasticity_domgroup[[i]])  %>% 
  mutate(estimate = estimate*attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:scale')+attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:center'),
         lower = lower*attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:scale')+attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:center'),
         upper = upper*attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:scale')+attr(scale(list_flow_plasticity_2[[i]]$perc_raw), 'scaled:center')
  )


list_dat[[i]] %>% 
  group_by(domgroup,housing) %>% 
  summarise(mean = mean(perc,na.rm = T),
            sd = sd(perc, na.rm = T),
            number = n()) %>% 
  mutate(sem = sd/sqrt(number-1)) -> domgroup_sum 
result_list_plasticity_domgroup[[i]] %>% summary -> sum



plasticity_plot_list[[i]] <- list_dat[[i]] %>% 
  ggplot(aes(housing, perc, group = subjectID)) +
  geom_line(aes(group=subjectID),color="grey80",size=0.8)+
  geom_point(color = "grey",fill = "grey80", shape=21,size=2, alpha = 0.7)+
  labs(x="",
       y="% of cells")+
  geom_ribbon(data = domgroup_sum,aes(x = housing, y = mean, ymin = mean - sem, ymax = mean + sem, fill = domgroup,group = domgroup), alpha = 0.3)+
  geom_line(data = domgroup_sum,aes(x = housing, y = mean, color = domgroup,group = domgroup),size = 1.5, alpha = 0.8)+
  geom_point(data = domgroup_sum,aes(x = housing, y = mean, color = domgroup,group = domgroup, fill = domgroup),size = 4, shape = 18, alpha =0.8)+
  scale_x_discrete(labels=c("PD09","GD14"))+
  scale_color_manual(values = c("purple4","#21908CFF","orange"))+
  scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
  ggtitle(glue('{names(list_dat)[[i]]} \n {round(sum$fixed[2,1],2)} [{round(sum$fixed[2,3],2)}, {round(sum$fixed[2,4],2)}]'))+
  theme_base(base_size = 9)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11),
        plot.background = element_blank())+
  annotate("pointrange",
         x = 2.2 ,
         y = cond_dat[[i]] %>% filter(housing=="group") %>% .$estimate, 
         ymin = cond_dat[[i]] %>% filter(housing=="group") %>% .$lower, 
         ymax = cond_dat[[i]] %>% filter(housing=="group") %>% .$upper,
         size = 0.7,
         alpha=0.9)+
  annotate("pointrange",
           x = 0.8, 
           y = cond_dat[[i]] %>% filter(housing=="pair") %>% .$estimate, 
           ymin = cond_dat[[i]] %>% filter(housing=="pair") %>% .$lower, 
           ymax = cond_dat[[i]] %>% filter(housing=="pair") %>% .$upper,
           size = 0.7,
           alpha=0.9) #+
  # annotate("text", x = 1.5, y = Inf-4, vjust = "inward", size = 4,
  #          label = glue('{round(sum$fixed[2,1],2)} [{round(sum$fixed[2,3],2)}, {round(sum$fixed[2,4],2)}]'))
}

# # for manuscript - b cell 
# list_dat[[1]] %>% 
#   group_by(housing, domgroup) %>% 
#   summarise(percs = mean(perc, na.rm =T)) %>% 
#   spread(housing, percs) %>% 
#   mutate(change = pair - group)



list_dat[[i]] %>% 
  ggplot(aes(housing, perc, group = subjectID)) +
  geom_line(aes(group=subjectID),color="grey80",size=0.8)+
  geom_point(color = "grey",fill = "grey80", shape=21,size=2.5, alpha = 0.6)+
  labs(x="",
       y="% of cells",
       color = "Social status",
       fill = "Social status")+
  geom_ribbon(data = domgroup_sum,aes(x = housing, y = mean, ymin = mean - sem, ymax = mean + sem, fill = domgroup,group = domgroup), alpha = 0.3)+
  geom_line(data = domgroup_sum,aes(x = housing, y = mean, color = domgroup,group = domgroup),size = 1, alpha = 0.7)+
  geom_point(data = domgroup_sum,aes(x = housing, y = mean, color = domgroup,group = domgroup, fill = domgroup),size = 6, shape = 18, alpha =0.8)+
  scale_x_discrete(labels=c("PD09","GD14"))+
  scale_color_manual(values = c("purple4","#21908CFF","orange"))+
  scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
  ggtitle(glue('{names(list_dat)[[i]]}'))+
  theme_base(base_size = 15)+
  theme(legend.position = "right",
        plot.title = element_text(size = 15),
        plot.background = element_blank()) -> for_legend

legend <- cowplot::get_legend(for_legend)
grid.newpage()
grid.draw(legend)



# plasticity_plot_list
# 
# main <- grid.arrange(plasticity_plot_list[[1]],plasticity_plot_list[[2]],plasticity_plot_list[[3]],
#                      plasticity_plot_list[[8]],plasticity_plot_list[[9]],plasticity_plot_list[[11]],
#                      plasticity_plot_list[[13]],legend, nrow = 2)
# 
# ggsave(main, filename = "results_figures/Figure4.png",width = 13,height = 10)



png(filename = glue("results_figures/supp_plasticity_new.png"),
    width = 25, height = 24, units = "cm", res = 600)


all <- grid.arrange(plasticity_plot_list[[1]],plasticity_plot_list[[3]],
                    plasticity_plot_list[[5]],plasticity_plot_list[[6]],plasticity_plot_list[[7]],
                    plasticity_plot_list[[8]],plasticity_plot_list[[9]],plasticity_plot_list[[10]],
                    plasticity_plot_list[[11]],plasticity_plot_list[[12]],
                     plasticity_plot_list[[13]],plasticity_plot_list[[14]],legend, ncol = 5)

invisible(dev.off())


do# ggsave(all, filename = "results_figures/supp_plasticity.png",width = 15,height = 14)
# 

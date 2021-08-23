
alldata %>%
  filter(!is.na(AGD)) %>%
  filter(AGD>10) %>% 
  ggplot(.,aes(AGD,ds))+
  geom_point()+
  facet_wrap(~km3_cluster)+
  geom_smooth(method = "lm")



AGD_ds <- brm(ds ~ AGD + (1|cohort), 
              data= alldata, file = "results_statRDS/AGD_ds.Rds")

AGD_ds # NO EFFECT  


AGD_ds1 <- brm(ds ~ AGD*km3_cluster + (1|cohort), 
              data= alldata, file = "results_statRDS/AGD_ds1.Rds")




alldata %>%
  filter(!is.na(AGD)) %>%
  filter(AGD>10) %>% 
  ggplot(.,aes(km3_cluster,AGD))+
  geom_point()+
  coord_flip()


alldata %>%
  filter(AGD>10) %>% 
  ggplot(.,aes(km3_cluster,AGD))+
  geom_point()+
  coord_flip()




bw_mod <- brm(yvar ~ ds + (1|cohort), 
              data= alldata %>% 
                mutate(yvar = scale(GD14_bw)),
              file = "results_statRDS/GD14_bw_ds.Rds")


bw_mod %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Bodyweight") -> plotdf_gd14_bw

# saveRDS(plotdf_gd14_bw,"results_figures/plotdf_gd14_bw.RDS")

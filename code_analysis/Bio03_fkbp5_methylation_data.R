

head(fkbp5)

fkbp5_df<- fkbp5 %>% 
  filter(timepoint =="GD14") %>% 
  mutate(cohort = LETTERS[cohort - 104]) %>% 
  mutate(subjectID = as.character(glue('{cohort}.{mouseID}'))) %>% 
  select(subjectID, CpG_1,CpG_2,CpG_3,Mean) %>% 
  left_join(.,alldata) %>% 
  left_join(.,resultsdf) %>% 
  filter(CpG_1>25) %>% 
  filter(CpG_2>25) %>% 
  filter(CpG_3>25) %>% 
  mutate(bwchange = GD14_bw-GD01_BW)  %>% 
  mutate(growth = bwchange/GD01_BW) %>% 
  mutate(growth_raw = growth) %>% 
  # mutate(CpG_1 = scale(CpG_1)) %>% 
  # mutate(CpG_2 = scale(CpG_3)) %>% 
  # mutate(CpG_3 = scale(CpG_1)) %>% 
  mutate(ds_rank = factor(ds_rank,
                          levels = c("1","2","3","4","5","6","7","8","9","10"))) %>%
  mutate(glicko_rank = factor(glicko_rank,
                          levels = c("1","2","3","4","5","6","7","8","9","10"))) %>% 
  mutate(growth = scale(growth))

 

fkbp5_df %>% 
  ggplot(.,aes(CpG_1))+
  geom_histogram()

fkbp5_df %>% 
  ggplot(.,aes(CpG_2))+
  geom_histogram()

fkbp5_df %>% 
  ggplot(.,aes(CpG_3))+
  geom_histogram()

fkbp5_df %>% 
  ggplot(.,aes(Mean))+
  geom_histogram()



fkbp5_df %>% 
  filter(CpG_3 <80)


fkbp5_df %>%
  gather(site,value,2:4) %>%
  ggplot(.,aes(site,value,group = subjectID,color=ds))+
  geom_line()+
  scale_color_viridis()


colnames(fkbp5_df)


# my_variable = 'glicko_rank'
plot_fkbp5 <- function(my_variable){
  temp <- fkbp5_df[,c(my_variable,'CpG_1','CpG_2','CpG_3')]
  colnames(temp)[1] <- 'my_variable'
  temp %>% 
    gather(key,value,2:4) %>% 
    filter(abs(value)>40) %>%
    ggplot(aes(my_variable,value))+
    geom_point(shape = 21, size = 1.5)+
    geom_point(shape = 21, size = 1.5,  fill = "grey", alpha = 0.5)+
    facet_wrap(~key, 
               labeller = labeller(key = c("CpG_1" ="CpG site 1",
                                          "CpG_2" = "CpG site 2",
                                          "CpG_3" = "CpG site 3")))+
    scale_y_continuous(breaks = c(40,60,80,100), limits = c(38,100))+
    labs(x = glue("{my_variable}"),
         y = "Fkbp5 DNA methylation (%)")+
    theme_bw(base_size = 6)
}

png(filename = glue("results_figures/Fkbp5_davidscore_biggerstripfont.png"),
    width = 10, height = 5, units = "cm", res = 600)
plot_fkbp5('ds')+xlab("David's score")+
  theme(  strip.text = element_text(size = rel(1.02)))
invisible(dev.off())

png(filename = glue("results_figures/Fkbp5_ds_rank.png"),
    width = 10, height = 5, units = "cm", res = 600)
plot_fkbp5('ds_rank')+xlab("Social rank")
invisible(dev.off())


png(filename = glue("results_figures/Fkbp5_cortpost.png"),
    width = 10, height = 5, units = "cm", res = 600)
plot_fkbp5('cort_post')+xlab("Corticosterone (GD14) (ng/ml)")
invisible(dev.off())


png(filename = glue("results_figures/Fkbp5_cortore.png"),
    width = 10, height = 5, units = "cm", res = 600)
plot_fkbp5('cort_pre')+xlab("Corticosterone (PD09) (ng/ml)")
invisible(dev.off())



plot_fkbp5('glicko_rank')
plot_fkbp5('cort_post')
plot_fkbp5('cort_pre')
plot_fkbp5('bwchange') # this is actually interesting
plot_fkbp5('growth') # of course, similar pattern as above 






plot_fkbp5_color <- function(my_variable){
  temp <- fkbp5_df[,c(my_variable,'CpG_1','CpG_2','CpG_3', "glicko_rank","ds")]
  colnames(temp)[1] <- 'my_variable'
  temp %>% 
    gather(key,value,2:4) %>% 
    # filter(abs(value)<3) %>% 
    ggplot(aes(my_variable,value,
               color = as.numeric(as.character(glicko_rank)),
               fill = as.numeric(as.character(glicko_rank))))+
    geom_point(alpha = 0.3,size = 2)+
    facet_wrap(~key, 
               labeller = labeller(key = c("CpG_1" ="CpG site 1",
                                           "CpG_2" = "CpG site 2",
                                           "CpG_3" = "CpG site 3")))+
    scale_color_distiller(palette = "Spectral")+
    scale_fill_distiller(palette = "Spectral")+
    # scale_color_gradient(high = "orange",low = "purple4")+
    # scale_fill_gradient(high = "orange",low = "purple4")+
    labs(x = glue("{my_variable}"),
         y = "Fkbp5 DNA methylation (%)",
         color = "Glicko rank",
         fill = "Glicko rank")+
    theme_bw(base_size = 6)+
    theme(legend.position = "top")
    
}

png(filename = glue("results_figures/Fkbp5_cortpre_col.png"),
    width = 10, height = 6, units = "cm", res = 600)
plot_fkbp5_color('cort_pre')+xlab("Corticosterone (PD09) (ng/ml)")
invisible(dev.off())


png(filename = glue("results_figures/Fkbp5_cortpost_col.png"),
    width = 10, height = 6, units = "cm", res = 600)
plot_fkbp5_color('cort_post')+xlab("Corticosterone (GD14) (ng/ml)")
invisible(dev.off())


plot_fkbp5('cort_post')

plot_fkbp5_color(my_variable = "cort_post") + xlab("Corticosterone (GD14) (ng/ml)")



plot_fkbp5('km3_cluster')
plot_fkbp5('ds_rank')
plot_fkbp5('glicko_rank')
plot_fkbp5('Win')
plot_fkbp5('Loss') #interesting 
plot_fkbp5('loss_prop') #interesting 

plot_fkbp5('ds')
plot_fkbp5('despotism')
plot_fkbp5('ind_cert')
plot_fkbp5('freeze_perc')
plot_fkbp5('spleen_weight_mg')
plot_fkbp5('pair_status')
plot_fkbp5('Lag')
plot_fkbp5('flee_perc')
plot_fkbp5('sub_pos_perc')
plot_fkbp5('bwchange') # this is actually interesting 
plot_fkbp5('growth') # of course, similar pattern as above 

# STAT ==============================================================================

# Null model 
cpg2_mod_null<-brm(CpG_2~1+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpg2_mod_null_95.Rds")


# GROWTH 
cpg1_mod_growth<-brm(CpG_1~growth+(1|cohort),
                     control = list(adapt_delta = 0.95),
                     data=fkbp5_df, file = "results_statRDS/cpg1_mod_growth.Rds")
cpg2_mod_growth<-brm(CpG_2~growth+(1|cohort),
                     control = list(adapt_delta = 0.95),
                     data=fkbp5_df, file = "results_statRDS/cpg2_mod_growth.Rds")
cpg3_mod_growth<-brm(CpG_3~growth+(1|cohort),
                     control = list(adapt_delta = 0.95),
                     data=fkbp5_df, file = "results_statRDS/cpg3_mod_growth.Rds")
cpgmean_mod_growth<-brm(Mean~growth+(1|cohort),
                        control = list(adapt_delta = 0.95),
                     data=fkbp5_df, file = "results_statRDS/cpgmean_mod_growth.Rds")


# growth and David's score 
cpg2_mod_growth_ds<-brm(CpG_2~growth+ds+(1|cohort),
                        control = list(adapt_delta = 0.95),
                     data=fkbp5_df, file = "results_statRDS/cpg2_mod_growth_ds.Rds")

# Loss - as cohort is a random effect no need to normalize 
cpg1_mod_loss<-brm(CpG_1~Loss+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpg1_mod_loss_95.Rds")


cpg2_mod_loss<-brm(CpG_2~Loss+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpg2_mod_loss_95.Rds")


cpg3_mod_loss<-brm(CpG_3~Loss+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpg3_mod_loss_95.Rds")


cpgmean_mod_loss<-brm(CpG_2~Loss+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpgmean_mod_loss_95.Rds")



# Loss and growth 
cpg2_mod_loss_growth<-brm(CpG_2~Loss+growth+(1|cohort),
                   control = list(adapt_delta = 0.95),
                   data=fkbp5_df, file = "results_statRDS/cpg2_mod_loss_growth_95x.Rds")

cpg2_mod_loss_growth_int<-brm(CpG_2~Loss*growth+(1|cohort),
                          control = list(adapt_delta = 0.95),
                          data=fkbp5_df, file = "results_statRDS/cpg2_mod_loss_growth_int_95.Rds")

cpg2_mod_loss_growth


# dominance certainty - no effect 
cpg1_mod_ind_cert<-brm(CpG_1~ind_cert+(1|cohort),
                       control = list(adapt_delta = 0.95),
                       data=fkbp5_df, file = "results_statRDS/cpg1_mod_ind_cert.Rds")
cpg2_mod_ind_cert<-brm(CpG_2~ind_cert+(1|cohort),
                       control = list(adapt_delta = 0.95),
                       data=fkbp5_df, file = "results_statRDS/cpg2_mod_ind_cert.Rds")
cpg3_mod_ind_cert<-brm(CpG_3~ind_cert+(1|cohort),
                       control = list(adapt_delta = 0.95),
                       data=fkbp5_df, file = "results_statRDS/cpg3_mod_ind_cert.Rds")
cpgmean_mod_ind_cert<-brm(Mean~ind_cert+(1|cohort),
                          control = list(adapt_delta = 0.95),
                          data=fkbp5_df, file = "results_statRDS/cpgmean_mod_ind_cert.Rds")

# Model comparison 
cpg2_mod_null <- add_criterion(cpg2_mod_null, "loo")
cpg2_mod_growth <- add_criterion(cpg2_mod_growth, "loo")
cpg2_mod_loss <- add_criterion(cpg2_mod_loss, "loo")
cpg2_mod_loss_growth <- add_criterion(cpg2_mod_loss_growth, "loo")


loo(cpg2_mod_null)
loo(cpg2_mod_growth)
loo(cpg2_mod_loss)
loo(cpg2_mod_loss_growth) #lowest loo value 

loo_compare(cpg2_mod_growth, cpg2_mod_loss_growth, criterion = "loo")

cpg2_mod_loss_growth # this model 
# plot(conditional_effects(cpg2_mod_loss_growth))



fkbp5_df %>% 
  ggplot(aes(Loss,growth,color = CpG_2, fill = CpG_2))+
  geom_point(size = 3, shape = 21, alpha =.5)+
  theme_bw()+
  scale_color_distiller(palette = 6)+
  scale_fill_distiller(palette = 6)



fkbp5_df %>% 
  ggplot(.,aes(Loss,CpG_2))+
  geom_point()
  
fkbp5_df %>% 
  filter(cohort!="C") %>% 
  ggplot(.,aes(loss_prop,CpG_1))+
  geom_point()+
  geom_smooth(method= "lm")



# David's score itselft 
cpg1_mod_ds<-brm(CpG_1~ds+(1|cohort),
                       data=fkbp5_df, file = "results_statRDS/cpg1_mod_ds.Rds")
cpg2_mod_ds<-brm(CpG_2~ds+(1|cohort),
                       data=fkbp5_df, file = "results_statRDS/cpg2_mod_ds.Rds")
cpg3_mod_ds<-brm(CpG_3~ds+(1|cohort),
                       data=fkbp5_df, file = "results_statRDS/cpg3_mod_ds.Rds")
# cpgmean_mod_ds<-brm(Mean~ds+(1|cohort),
#                           data=fkbp5_df, file = "results_statRDS/cpgmean_mod_ds.Rds")



cpg1_mod_ds %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Fkbp5 CpG1") -> x

cpg2_mod_ds %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Fkbp5 CpG2") -> y 

cpg3_mod_ds %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Fkbp5 CpG3") -> z 

rbind(x,y,z) -> plotdf_gd14_fkbp5
# saveRDS(plotdf_gd14_fkbp5,"results_statRDS/plotdf_gd14_fkbp5.RDS")

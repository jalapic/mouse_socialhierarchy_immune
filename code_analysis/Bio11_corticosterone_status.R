alldata %>% 
  select(batch,pairID,cohort,subjectID,mouseID,cort_pre,cort_post,
         domgroup,ds_rank,ds,win_prop,loss_prop,pair_status,despotism) %>% 
  mutate(domgroup = factor(domgroup, levels = c("Alpha","Subdominant","Subordinate"),ordered = T)) %>% 
  mutate(despotism_binary = ifelse(despotism >=0.5,"Despotism: High","Despotism: Low")) -> cortdf

# saveRDS(cortdf, "results_statRDS/cortdf.RDS")

cortdf %>% 
  mutate(buddyID = paste(batch, mouseID, sep = "-")) %>% 
  ggplot(aes(pair_status, cort_pre, group = buddyID))+
  geom_point()+
  geom_line()


cortdf %>% 
  ggplot(aes(ds, cort_post))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~cohort) # makes me wonder about despotism 


cortdf %>% 
  ggplot(aes(cort_pre,cort_post))+
  geom_point() # this is actually somewhat interesting 

cortdf %>% 
  gather(key,value, c('cort_pre','cort_post')) %>% 
  mutate(key = ifelse(key == "cort_post", "Post","Pre")) %>% 
  mutate(key = factor(key, levels = c('Pre','Post'))) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", "Dominant","Subordinate")) %>% 
  ggplot(aes(key,value,group = subjectID, color = pair_status, fill = pair_status))+
  geom_line(alpha = 0.5)+
  geom_point(shape = 21, alpha = 0.5,)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~domgroup)+
  theme(legend.position = "top")+
  labs(x = "", y = "Corticosterone level (ng/ml)",
    color = 'Pair housing social status',
       fill = 'Pair housing social status')
# so run mixed effects model - tiempoint*domgroup 

cortdf %>% 
  ggplot(aes(domgroup,cort_post))+
  geom_boxjitter(outlier.color = NA)

all_behavior %>% 
  select(cohort, despotism) %>% 
  unique()


cortdf %>% 
  ggplot(aes(domgroup,cort_post, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none") +
  labs(x = "Social status", y = "Corticosterone level (ng/ml)")+
  facet_wrap(~despotism_binary) # awesome, pretty much replicating the hormone paper


cortdf %>%
  mutate(despotism_binary = ifelse(despotism >=0.5,"Despotism: High","Despotism: Low")) %>%
  ggplot(aes(ds,cort_post))+
  geom_point()+
  facet_wrap(~despotism_binary) # domgroup should be better

cortdf %>% 
  filter(domgroup == "Alpha") %>% 
  ggplot(aes(despotism, cort_post))+
  labs(title = "Only alphas, text = cohort ID")+
  geom_text(aes(label = cohort))  # interesting that it's not negatively correlated 



cortdf %>% 
  ggplot(aes(domgroup,cort_post, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none") +
  labs(x = "Social status", y = "Corticosterone level (ng/ml)")+
  facet_wrap(~pair_status) 


## now stats ============================================================================
# pair housing data 
cort_pre <- brm(cort_pre ~ pair_status + (1|pairID),
                    data = cortdf %>% 
                      mutate(pairID = paste(batch,mouseID,sep = "_")),
                    prior = set_prior('normal(0, 30)'),
                    file = "results_statRDS/cort_pre.Rds")

cort_pre
marginal_effects(cort_pre)



## Group housing corticosterone level - data generating process 
cort_ds_mod0 <- brm(cort_post ~ 1 + (1|cohort),
                    data = cortdf,
                    prior = set_prior('normal(0, 30)'),
                    iter = 1000, chains = 4, cores = 4,
                    file = "results_statRDS/cort_ds_mod0.Rds")


cort_ds_mod1 <- brm(cort_post ~ ds + (1|cohort),
                       data = cortdf,
                       file = "results_statRDS/cort_ds_mod1.Rds")


cort_ds_mod2 <- brm(cort_post ~ ds*despotism_binary + (1|cohort),
                  data = cortdf,
                  file = "results_statRDS/cort_ds_mod2.Rds")

cort_ds_mod3 <- brm(cort_post ~ mo(domgroup) + (1|cohort),
                    data = cortdf,
                    file = "results_statRDS/cort_ds_mod3.Rds")

# estimates  
summary(cort_ds_mod0)
summary(cort_ds_mod1)
summary(cort_ds_mod2)
summary(cort_ds_mod3)

# model diagnostics 
plot(cort_ds_mod0)
plot(cort_ds_mod1)
plot(cort_ds_mod2)



# compare three models 

cort_loo_null <- add_criterion(cort_ds_mod0, "loo")
cort_loo_ds <- add_criterion(cort_ds_mod1, "loo")
cort_loo_despotism <- add_criterion(cort_ds_mod2, "loo")

loo_compare(cort_loo_null, cort_loo_ds, cort_loo_despotism, criterion = "loo")

# loo(cort_loo_null)
# loo(cort_loo_ds)
# loo(cort_loo_despotism)


cort_mod_ds <- brm(yvar ~ ds + (1|cohort),
                    data = cortdf %>% 
                     mutate(yvar = scale(cort_post)),
                    file = "results_statRDS/cort_mod_ds_scaled.Rds")

cort_mod_ds %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Corticosterone") -> plotdf_gd14_cort

# saveRDS(plotdf_gd14_cort,"results_figures/plotdf_gd14_cort.RDS")



# low and high despotism comparison separately 
summary(lme4::lmer(cort_post ~ domgroup + (1|cohort),
           data = cortdf %>% 
             mutate(domgroup = factor(domgroup, levels = c("Subordinate","Alpha","Subdominant"))) %>% 
             filter(despotism_binary =="Despotism: High")))

cort_domgroup_despot_high <- brm(cort_post ~ 0 + domgroup + (1|cohort),
                    data = cortdf %>% 
                      filter(despotism_binary =="Despotism: High"),
                    file = "results_statRDS/cort_domgroup_high_despot.Rds")
cort_domgroup_despot_high
conditional_effects(cort_domgroup_despot_high)

cort_domgroup_despot_low <- brm(cort_post ~ 0 + domgroup + (1|cohort),
                                 data = cortdf %>%  
                                   filter(despotism_binary !="Despotism: High"),
                                 file = "results_statRDS/cort_domgroup_low_despot.Rds")

cort_domgroup_despot_low


get_variables(cort_domgroup_despot_high)

get_variables(cort_domgroup_despot_low)

cort_domgroup_despot_high %>% 
  conditional_effects() %>% .$domgroup
  as.data.frame()

cort_domgroup_despot_low %>% 
  spread_draws(b_domgroupAlpha, b_domgroupSubdominant, b_domgroupSubordinate) %>% 
  mutate(alpha_subdom = b_domgroupAlpha - b_domgroupSubdominant,
         alpha_sub = b_domgroupAlpha - b_domgroupSubordinate,
         subdom_sub = b_domgroupSubdominant - b_domgroupSubordinate) %>% 
  median_qi() %>% 
  t
  

cort_domgroup_despot_high %>% 
  spread_draws(b_domgroupAlpha, b_domgroupSubdominant, b_domgroupSubordinate) %>% 
  mutate(alpha_subdom = b_domgroupAlpha - b_domgroupSubdominant,
         alpha_sub = b_domgroupAlpha - b_domgroupSubordinate,
         subdom_sub = b_domgroupSubdominant - b_domgroupSubordinate) %>% 
  select(alpha_subdom, alpha_sub, subdom_sub) %>% 
  median_qi() %>% t


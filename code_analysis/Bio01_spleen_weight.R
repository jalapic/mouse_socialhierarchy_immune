
colnames(alldata)


# David's score
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(ds,spleen_bw_ratio))+
  labs(x="David's score",
       y="Spleen weight/body weight(mg/g)")+
  geom_point(size = 2, shape = 21)+
  theme_bw()+
  geom_vline(xintercept = 0, color = "red", linetype = "dashed")


# domgroup
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(domgroup,spleen_bw_ratio))+
  labs(x="Social status group",
       y="Spleen weight/body weight(mg/g)")+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 jitter.height = 0.01, jitter.width = 0.1, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  theme_base()



# behavior cluster 
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(km3_cluster,spleen_bw_ratio))+
  labs(x="Cluster",
       y="Spleen weight/body weight(mg/g)")+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 jitter.height = 0.01, jitter.width = 0.1, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  theme_base()


# aggression received 
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(loss_prop,spleen_bw_ratio))+
  geom_point()+
  theme_bw()


# aggression received 
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(Mono_GD14,spleen_bw_ratio))+
  geom_point()+
  theme_bw()+
  facet_wrap(~km3_cluster)


# aggression given - meh 
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(Win,spleen_bw_ratio))+
  geom_point()+
  geom_smooth(method = "lm")


alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(pair_status,spleen_bw_ratio))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 jitter.height = 0.01, jitter.width = 0.1, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  theme_base()


# defensive behaviors - all meh 
alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(freeze_perc,spleen_bw_ratio))+
  geom_point()+
  geom_smooth(method = "lm")


alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(sub_pos_perc,spleen_bw_ratio))+
  geom_point()+
  geom_smooth(method = "lm")


alldata %>%
  filter(!is.na(spleen_weight_mg)) %>%
  filter(!is.na(GD14_bw)) %>%
  mutate(spleen_bw_ratio=spleen_weight_mg/GD14_bw) %>%
  ggplot(.,aes(flee_perc,spleen_bw_ratio))+
  geom_point()+
  geom_smooth(method = "lm")




#spleen - brms linear fit ====================================================================
sp_mod2<-brm(spleen_weight_mg/GD14_bw~ds+(1|cohort),data=alldata,file = "results_statRDS/spleen_davidscore.Rds")

summary(sp_mod2)

dat<-conditional_effects(sp_mod2)$ds
dat2<-dat %>% mutate(xvar=ds,yvar=estimate__,lower=lower__,upper=upper__) %>% select(xvar,yvar,lower,upper) %>% as.data.frame()

sp<-ggplot()+
  geom_point(data=alldata,
             aes(ds,spleen_weight_mg/GD14_bw),
             alpha = 0.3,shape=21,size=3.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="David's score",
       y="Spleen weight/Body weight ratio")+
  theme_bw()+
  theme(legend.position = c(0.85,0.85))+
  annotate("text", x = 7, y = 4.3, label = "Posterior median [95% CI] = -0.01 [-0.02, -0.00]", size = 3)

print(sp)

ggsave(sp, filename = "results_figures/spleen_ds.png",width = 6,height = 5)



sp_mod_ds<-brm(yvar~ds+(1|cohort),
               data=alldata %>% 
                 mutate(yvar = scale(spleen_weight_mg/GD14_bw)),
               file = "results_statRDS/spleen_davidscore_scaled.Rds")


sp_mod_ds %>% 
  spread_draws(b_ds) %>% #change
  median_qi(.width = c(.95, .66)) %>% 
  mutate(cell_type = "Spleen weight") -> plotdf_gd14_spleen

# saveRDS(plotdf_gd14_spleen,"results_figures/plotdf_gd14_spleen.RDS")




sp_mod_rank<-brm(spleen_weight_mg/GD14_bw~mo(glicko_rank)+(1|cohort),data=alldata,file = "results_statRDS/spleen_moglick.Rds")

summary(sp_mod_rank)
conditional_effects(sp_mod_rank)


# despotism ==================================================================================

library(mice)

md.pattern(all_behavior) # impute missing bodyweight value for one alpha male 
all_behaviorx <- mice(all_behavior,m=5,maxit=50,meth='pmm',seed=500, verbose = F)



alphadata <- alldata %>% 
  filter(glicko_rank == 1) %>% 
  select(subjectID, despotism, GD14_bw, spleen_weight_mg)

alphadata [2:3,3] <- all_behaviorx$imp$GD14_bw[,1]
  
sp_mod3<-brm(spleen_weight_mg/GD14_bw~ despotism, 
             data = alphadata, 
             file = "results_statRDS/sp_mod_despotism.Rds")
summary(sp_mod3)

dat<-conditional_effects(sp_mod3)$despotism
dat2<-dat %>% mutate(xvar=despotism,yvar=estimate__,lower=lower__,upper=upper__) %>% select(xvar,yvar,lower,upper) %>% as.data.frame()



spx<-ggplot()+
  geom_point(data=alphadata,
             aes(despotism,spleen_weight_mg/GD14_bw),
             alpha = 0.7,shape=21,size=3.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="Despotism",
       y="Spleen weight/Body weight ratio",
       color = "Cluster",
       fill = "Cluster")+
  theme_bw()+
  theme(legend.position = c(0.85,0.85))+
  annotate("text", x = 0.56, y = 2.48, label = "Posterior median [95% CI] = 2.71 [0.94, 4.55]", size = 3)

print(spx)

ggsave(spx, filename = "results_figures/spleen_despot.png",width = 6,height = 5)

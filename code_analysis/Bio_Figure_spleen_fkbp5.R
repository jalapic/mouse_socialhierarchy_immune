
dat<-conditional_effects(sp_mod2)$ds
dat2<-dat %>% mutate(xvar=ds,yvar=estimate__,lower=lower__,upper=upper__) %>% select(xvar,yvar,lower,upper) %>% as.data.frame()

sp<-ggplot()+
  geom_point(data=alldata,
             aes(ds,spleen_weight_mg/GD14_bw),
             shape = 21, size = 1.5)+
  geom_point(data=alldata,
             aes(ds,spleen_weight_mg/GD14_bw),
             shape = 21, size = 1.5,  fill = "grey", alpha = 0.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="David's score",
       y="Spleen weight/Body weight ratio",
       color = "Social status group",
       fill = "Social status group")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position ="none")+
  annotate("text", x = 4.55, y = 6.5, label = "Posterior median [95% CI] = -0.10 [-0.16, -0.03]", 
           size = 2, color = "red")+
  theme_bw(base_size = 6)


png(filename = glue("results_figures/spleen_weight_davidscore.png"),
    width = 6, height = 5, units = "cm", res = 600)
sp
invisible(dev.off())




dat<-conditional_effects(sp_mod3)$despotism
dat2<-dat %>% mutate(xvar=despotism,yvar=estimate__,lower=lower__,upper=upper__) %>% select(xvar,yvar,lower,upper) %>% as.data.frame()


spx<-ggplot()+
  geom_point(data=alphadata,
             aes(despotism,spleen_weight_mg/GD14_bw),
             shape = 21, size = 1.5)+
  geom_point(data=alphadata,
             aes(despotism,spleen_weight_mg/GD14_bw),
             shape = 21, size = 1.5,  fill = "grey", alpha = 0.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="Despotism",
       y="Spleen weight/Body weight ratio",
       color = "Social status group",
       fill = "Social status group")+
  theme_bw(base_size = 6)+
  annotate("text", x = 0.5, y = 2.08, 
           label = "Posterior median [95% CI] = 2.71 [0.94, 4.55]", size = 2, color = 'red')


png(filename = glue("results_figures/spleen_weight_despotism.png"),
    width = 6, height = 5, units = "cm", res = 600)
spx
invisible(dev.off())




# CpG_2 ~ growth

dat <- conditional_effects(cpg2_mod_loss_growth)$growth
x = 2
mybackscale <- function (x){
  t <- x * attr(scale(fkbp5_df$growth_raw), 'scaled:scale') + attr(scale(fkbp5_df$growth_raw), 'scaled:center') 
  return (t)}

dat2 <- dat %>% 
  mutate(xvar=growth,yvar=estimate__,lower=lower__,upper=upper__) %>% 
  select(xvar,yvar,lower,upper) %>% 
  as.data.frame() 
mybackscale(dat2$xvar) -> dat2$xvar


cpg2_growth<-ggplot()+
  geom_point(data=fkbp5_df,
             aes(x = growth_raw*100, y = CpG_2,color = domgroup, fill = domgroup),
             alpha = 0.3,shape=21,size=3.5)+
  geom_ribbon(data=dat2,aes(xvar*100,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar*100,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="Body weight change over group housing period(%)",
       y="fkbp5 DNA methylation (CpG2 site)",
       color = "Social status group",
       fill = "Social status group")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  annotate("text", x = 1, y = 85, label = "2.40 [0.17, 4.70]", size = 4)



# CpG_2 ~ Loss 
dat <- conditional_effects(cpg2_mod_loss_growth)$Loss

dat2 <- dat %>% 
  mutate(xvar=Loss,yvar=estimate__,lower=lower__,upper=upper__) %>% 
  select(xvar,yvar,lower,upper) %>% 
  as.data.frame() 


cpg2_Loss<-ggplot()+
  geom_point(data=fkbp5_df,
             aes(x = Loss, y = CpG_2,color = domgroup, fill = domgroup),
             alpha = 0.3,shape=21,size=3.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="Number of losses over group housing period",
       y="fkbp5 DNA methylation (CpG2 site)",
       color = "Social status group",
       fill = "Social status group")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  annotate("text", x = 220, y = 85, label = "-0.04 [-0.08, -0.01]", size = 4)




#CpG1 ~ Loss 
dat <- conditional_effects(cpg1_mod_loss)$Loss

dat2 <- dat %>% 
  mutate(xvar=Loss,yvar=estimate__,lower=lower__,upper=upper__) %>% 
  select(xvar,yvar,lower,upper) %>% 
  as.data.frame() 


cpg1_Loss<-ggplot()+
  geom_point(data=fkbp5_df,
             aes(x = Loss, y = CpG_1,color = domgroup, fill = domgroup),
             alpha = 0.3,shape=21,size=3.5)+
  geom_ribbon(data=dat2,aes(xvar,ymax=upper,ymin=lower),alpha=0.1,fill="red")+
  geom_line(data=dat2,aes(xvar,yvar),size=1.5,alpha=0.5,color="red")+
  labs(x="Number of losses over group housing period",
       y="fkbp5 DNA methylation (CpG1 site)",
       color = "Cluster",
       fill = "Cluster")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  annotate("text", x = 220, y = 45, label = "-0.04 [-0.07, -0.01]", size = 4)



legend <- cowplot::get_legend(sp+
                                theme(legend.position = c(0.5,0.5)))

figure2 <- plot_grid(sp, spx, legend, cpg1_Loss, cpg2_Loss, cpg2_growth, labels = c("A","B","","C","D","E"), nrow = 2)


# Let's print 

theme_set(theme_bw(base_size = 10))

png(filename = "results_figures/flow_pairstatus.png",
    width = 18, height = 12, units = "cm", res = 600)
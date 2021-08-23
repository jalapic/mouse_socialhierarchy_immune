df.behav_freq %>%
  select(-mouseID, -cohort) %>%
  left_join(.,all_behavior %>% 
              select(domgroup,ds,subjectID,km3_cluster))  -> cluster_exp


df.behav_freq %>%
  left_join(hour) %>%
  select(-mouseID, -cohort) %>%
  mutate_if(is.numeric,~.x/totalobs) %>%
  select(-totalobs) %>%
  left_join(.,all_behavior %>% 
              select(domgroup,ds,subjectID,km3_cluster))  -> cluster_exp_hourly


myplot <- cluster_exp %>% 
  mutate(domgroup = factor(domgroup, levels = c("Subordinate","Subdominant","Alpha"))) %>% 
  group_by(domgroup, km3_cluster) %>%
  tally() %>% 
  ggplot(aes(km3_cluster,domgroup, label = as.character(n)))+
  geom_tile(aes( fill = n),alpha = 0.8)+
  geom_text(size = 14)+
  scale_fill_distiller(palette = "Spectral")+
  newggtheme+
  labs(x = "K-means clusters",
       y = "Social status group \n(based on David's score)")
ggsave(myplot, filename = "results_figures/contingency.png",width = 6,height = 5)

colnames(cluster_exp_hourly)

cluster_exp_hourly %>% 
  mutate(hourly_wins = lunge,fight,chase,mount) %>% 
  mutate(hourly_losses = preem_flee, flee, sub_pos,freeze) %>% 
  ggplot(aes(hourly_wins, hourly_losses,color = km3_cluster,fill = km3_cluster))+
  geom_point(size = 3.5,alpha = 0.35,shape = 21)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  labs(x = "Freq. of aggression given (per hour)",
       y = "Freq. of aggression received (per hour)",
       color = "Cluster",
       fill = "Cluster")+
  theme(legend.position = c(0.8,0.7)) -> myplot
ggsave(myplot, filename = "results_figures/cluster_agg_given_received.png",width = 6,height = 5)


cluster_exp_hourly %>% 
  ggplot(aes(ds,color = km3_cluster,fill = km3_cluster))+
  geom_histogram(alpha = 0.3)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  labs(x = "David's score",
       y = "Number of individuals",
       color = "Cluster",
       fill = "Cluster")+
  theme(legend.position = c(0.85,0.85)) -> myplot

 ggsave(myplot, filename = "results_figures/cluster_ds_histo.png",width = 6,height = 5)


bhv_df <- data.frame(
  behavior =  c("fight","mount","chase","lunge","flee","freeze","sub_pos","preem_flee","sniff","allogroom"),
  behaviorx = c("Fight","Mount","Chase","Lunge","Flee","Freeze","Subordinate posture","Preemtive flee","Sniff","Allogroom"))


cluster_exp_hourly %>% 
  select(-none,-walk) %>% 
  gather(behavior,value,1:10) %>% 
  left_join(bhv_df) %>% 
  mutate(behaviorx = factor(behaviorx,
                           levels = c("Fight","Mount","Chase","Lunge",
                                      "Flee","Freeze","Subordinate posture",
                                      "Preemtive flee","Sniff","Allogroom"))) %>%
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(km3_cluster,value,color = km3_cluster, fill = km3_cluster))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, 
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.09, errorbar.draw = TRUE,
                 position = position_dodge(0.85))+
  facet_wrap(~behaviorx, scales = "free_y", ncol = 5)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "Cluster",
       y = "Frequency (per hour)")+
  theme_bw()+
  theme(legend.position = "none") -> myplot
ggsave(myplot, filename = "results_figures/cluster_each_behavior.png",width = 12,height = 5)


#network measures - exploratory

indiv_network_df %>% 
  left_join(all_behavior %>% select(subjectID,ds_rank)) %>% 
  gather(net,value,1:7) %>% 
  ggplot(aes(ds_rank,value))+
  facet_wrap(~net, scales =  "free_y", ncol = 4)+
  geom_point()+
  geom_line(aes(color = cohort))+
  theme_bw()+
  scale_x_continuous(breaks = c(1:10))+
  scale_color_viridis(discrete = T) -> net_plot

ggsave(net_plot,filename = "results_figures/net_plot.png",height = 5, width = 11)


# Want to get some sort of 'loading score in k-means clustering' --> called feature detection 

                    
df <- df.behav_freq %>%
  select(-mouseID, -cohort) %>% 
  column_to_rownames(var = "subjectID") 


dis = dist(df)
sil = silhouette(km.res3$cluster, dis)
summary(sil)

png(filename = "results_figures/silhoutte_plot.png", 
    width = 800, height = 800, res = 120)

plot(sil)
dev.off()


set.seed(33)
res <- kcca(df,k=3)

FeatureImp_res <- FeatureImpCluster(res,as.data.table(df))


png(filename = "results_figures/kmeans_misclassification_rate.png", res = 100)
plot(FeatureImp_res)
dev.off()


png(filename = "results_figures/kmeans_barplot.png", 
    width = 400, height = 800, res = 120)
barplot(res)
dev.off()


# bwplot(res,df)

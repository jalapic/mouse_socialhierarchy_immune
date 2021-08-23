
colnames(all_behavior)

all_behavior %>% 
  #  filter(glicko_rank!=1) %>% 
  ggplot(.,aes(Rating,freeze_perc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Glicko rating",
       y="% of freeze")+
  theme_classic(base_size = 20)

all_behavior %>% 
  filter(glicko_rank!=1) %>% 
  ggplot(.,aes(Rating,flee_perc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Glicko rating",
       y="% of flee")+
  theme_classic(base_size = 20)

all_behavior %>% 
  filter(glicko_rank!=1) %>% 
  ggplot(.,aes(Rating,sub_pos_perc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Glicko rating",
       y="% of subordinate posture")+
  theme_classic(base_size = 20)



all_behavior %>% 
  filter(glicko_rank!=1) %>% 
  ggplot(.,aes(domgroup,sub_pos_perc))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x="Glicko rating",
       y="% of subordinate posture")+
  theme_classic(base_size = 20)



# multinomial - predictor day and aggresor rank
# 
# 
# multinom.mod1<-brm(response~domgroup+(1|cohort)+(1|subjectID),
#                    data=response_each_rank,
#                    family = categorical(),
#                    file = "results_statRDS/multinom.mod1.RDS",
#                    chains = 3,iter = 1000)
# 
multinom_mod2<-brm(response~mo(glicko_rank)+(1|cohort)+(1|subjectID),
                   data=response_each_rank,
                   family = categorical(link = "logit"),
                   file = "results_statRDS/multinom_mod2.RDS",
                   chains = 3,iter = 1000)


response_each_rankx<-response_each_rank %>% 
  mutate(response=factor(response,levels=c("Sub","Freeze","Flee")))

multinom_mod2x<-brm(response~mo(glicko_rank)+(1|cohort)+(1|subjectID),
                   data=response_each_rankx,
                   family = categorical(link = "logit"),
                   file = "results_statRDS/multinom_mod2x.RDS",
                   chains = 3,iter = 1000)



multinom_mod2xx<-brm(response~mo(glicko_rank)+(1|cohort)+(1|subjectID),
                     data=response_each_rankx,
                     family = categorical(link = "logit"),
                     file = "results_statRDS/multinom_mod2xx",
                     chains = 3,iter = 1000,
                     control = list(adapt_delta=0.95))



# Behavioral response stacked bar==============================

df<-all_behavior %>% 
  group_by(glicko_rank) %>% 
  summarise(Flee=mean(flee_perc,na.rm = T),
            `Subordinate Posture`=mean(sub_pos_perc,na.rm = T),
            Freeze=mean(freeze_perc,na.rm = T)) %>% 
  gather(Response,Proportion,2:4) %>% 
  mutate(glicko_rank=as.numeric(glicko_rank))


fig2 <- df %>% 
  ggplot(.,aes(glicko_rank,Proportion,fill=Response))+
  geom_bar(stat = "identity",alpha=0.78,position = "fill")+
  scale_fill_viridis(discrete=T)+
  newggtheme_with_legends+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  labs(x="Social rank",
       y="Proportion of responses to aggressive interactions")+
  theme(axis.title.y = element_text(size=15,vjust=2),
        axis.title.x = element_text(size=15))

print(fig2)


all_behavior %>% 
  group_by(glicko_rank,cohort) %>% 
  summarise(Flee=mean(flee_perc,na.rm = T),
            `Subordinate Posture`=mean(sub_pos_perc,na.rm = T),
            Freeze=mean(freeze_perc,na.rm = T)) %>%
  ungroup %>% 
  gather(Response,Proportion,3:5) %>% 
  mutate(glicko_rank=as.numeric(glicko_rank)) %>% 
  ggplot(.,aes(glicko_rank,Proportion,fill=Response))+
  geom_bar(stat = "identity",alpha=0.78,position = "fill")+
  scale_fill_viridis(discrete=T)+
  newggtheme_with_legends+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  labs(x="Social rank",
       y="Proportion of responses to aggressive interactions")+
  theme(axis.title.y = element_text(size=15,vjust=2),
        axis.title.x = element_text(size=15))+
  facet_wrap(~cohort)





all_behavior %>% colnames()
all_behavior %>% select(cohort,Loss,glicko_rank) %>% 
  group_by(cohort) %>% 
  mutate(loss_prop=Loss/sum(Loss)) %>% 
  ungroup() %>% 
  ggplot(.,aes(as.numeric(glicko_rank),loss_prop))+
  geom_point(size=5,shape=21)+
  newggtheme+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  labs(x="Social rank",y="Received aggression (%)")


all_behavior %>% select(cohort,Win,glicko_rank) %>% 
  group_by(cohort) %>% 
  mutate(win_prop=Win/sum(Win)) %>% 
  ungroup() %>% 
  ggplot(.,aes(as.numeric(glicko_rank),win_prop))+
  geom_point(size=5,shape=21)+
  newggtheme+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  labs(x="Social rank",y="Given aggression (%)")



p1 <- all_behavior %>% select(cohort,Loss,ds) %>% 
  group_by(cohort) %>% 
  mutate(loss_prop=Loss/sum(Loss)) %>% 
  ungroup() %>% 
  ggplot(.,aes(ds,loss_prop))+
  geom_point(size=3,shape=21)+
  newggtheme+
  labs(x="Normalized David's score",y="Received aggression (%)")


p2 <- all_behavior %>% select(cohort,Win,ds) %>% 
  group_by(cohort) %>% 
  mutate(win_prop=Win/sum(Win)) %>% 
  ungroup() %>% 
  ggplot(.,aes(ds,win_prop))+
  geom_point(size=3,shape=21)+
  newggtheme+
  labs(x="Normalized David's score",y="Given aggression (%)")

ds_agg <- grid.arrange(p2,p1, ncol = 2)
ggsave(ds_agg,file = "results_figures/ds_agg.png",height = 4, width = 8)


multinom_mod3<-brm(response~ds+(1|cohort)+(1|subjectID),
                   data=response_each_rank,
                   family = categorical(link = "logit"),
                   file = "results_statRDS/multinom_mod3.RDS",
                   chains = 3,iter = 1000)


response_each_rankx<-response_each_rank %>% 
  mutate(response=factor(response,levels=c("Sub","Freeze","Flee")))

multinom_mod3x<-brm(response~ds+(1|cohort)+(1|subjectID),
                    data=response_each_rankx,
                    family = categorical(link = "logit"),
                    file = "results_statRDS/multinom_mod3x.RDS",
                    chains = 3,iter = 1000)



multinom_mod3xx<-brm(response~ds+(1|cohort)+(1|subjectID),
                     data=response_each_rankx,
                     family = categorical(link = "logit"),
                     file = "results_statRDS/multinom_mod3xx",
                     chains = 3,iter = 1000,
                     control = list(adapt_delta=0.95))


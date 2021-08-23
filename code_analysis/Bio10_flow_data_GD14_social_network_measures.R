colnames(alldata)

flow_GD14 <- flow_all %>% 
  filter(housing=="group") %>% 
  filter(cohort!="Pair-comparison") %>% 
  left_join(alldata %>% select(subjectID, ds, inet.betw))


#
flow_GD14 %>% 
  ggplot(aes(inet.betw,perc))+
  geom_point()+
  facet_wrap(~ cell_type, scales = "free_y")


hist(flow_GD14$inet.betw)

#stat - won't be necessary, not going to include in the paper 



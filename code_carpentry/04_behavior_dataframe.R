# allbehav<-df.group %>%
#   bind_rows() %>%
#   mutate(behav=paste(`Animal 1 - Behavior`, `Animal 2 - Behavior`, sep="---"))
# 
# behavdf<- allbehav %>%
#   select(behav) %>%
#   unique()


## be careful not to overwrite the previous one - handcoded!
# behavdf_old<-read_csv("data_clean/behavior_score.csv")[,1:2]
# behavdf
# behavdf2 <- left_join(behavdf,behavdf_old)
# write.csv(behavdf2,"data_clean/behavior_scored.csv",row.names = F)

behavdf<-read_csv("data_clean/behavior_scored.csv")[,1:2]


df.group <-df.group %>% 
  map(~mutate(.,behav=paste(`Animal 1 - Behavior`, `Animal 2 - Behavior`, sep="---"))) %>% 
  map(~left_join(.,behavdf))

df.groupx<-lapply(df.group,get_df_new_ethogram)

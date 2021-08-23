## subordinate behavior pattern 

## How many times did the animal initiate social interaction? 
df.init <- df.group %>% 
  map(initiation) %>% 
  map2_df(.,names(.),~mutate(.x,cohort=.y)) %>% 
  mutate(subjectID=paste(cohort,mouseID,sep="."))


df.init



## Each behavior - frequency

l.each_behav <- df.group %>% map(each_behavior)

l.behav_freq <- l.each_behav %>% map(behavior_freq)

l.behav_perc <- l.behav_freq %>% map(behavior_perc)

df.behav_freq <- l.behav_freq %>% cohort_rbindlist() 
df.behav_perc <- l.behav_freq %>% map(behavior_perc) %>% 
  cohort_rbindlist()


l.sub_behav_perc <- l.behav_freq %>% map(sub_behavior_perc)

df.sub_behav_freq <- df.behav_freq %>% select(cohort,mouseID,preem_flee,flee,sub_pos,freeze)

df.sub_behav_perc <- l.sub_behav_perc %>% 
  cohort_rbindlist() %>% 
  mutate(mouseID=as.character(mouseID))






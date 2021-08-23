#01. observation descriptives
all_behav <- df.groupx %>%
  map2_df(.,names(.),~mutate(.x,cohort=.y)) %>% 
  select(cohort,day,hour,uniqueobs)

# number of days of observation
day<-all_behav %>%
  group_by(cohort) %>%
  dplyr::summarize(maxday = max(day)) %>%
  ungroup() %>%
  .$maxday %>% as.numeric() 
day
day%>%summary()

# number of observation hour
hour<-all_behav %>%
  group_by(cohort) %>%
  dplyr::summarize(totalobs = length(unique(uniqueobs)))
as.data.frame(hour)
hour%>%
  ungroup() %>%
  .$totalobs %>%
  sum

hour%>%
  summary()

#total aggressive contests
all_behav %>% 
  nrow()

#total hour PER day
all_behav %>%
  group_by(cohort) %>% 
  top_n(1,uniqueobs) %>% 
  as.data.frame() %>%
  summarise(ave_total_obs_hour=mean(uniqueobs/day))



df.groupx %>% 
  map(~select(.,Timestamp,winner,loser)) %>% 
  map(get_uncertainty) -> res


indivs <- lapply(res, function(x) x$indivs)

ind_cert_df<-indivs %>% 
  map2_df(.,names(.),~mutate(.x,cohort=.y)) %>% 
  mutate(subjectID=paste(cohort,ID,sep="."),
         ind_cert=Mean) %>% 
  select(cohort,subjectID,ind_cert) 



df.groupxx <- df.groupx %>% 
  map(~ select(., winner,loser,result)) %>% 
  map(~ as.data.frame(.))



indiv_network_df <- df.groupxx %>% 
  map(get_indiv_network_measures) %>% 
  map2_df(.,names(.),~mutate(.x,cohort = .y)) %>% 
  mutate(subjectID = paste(cohort,mouseID,sep = ".")) %>% 
  select(-mouseID)


str(indiv_network_df)


  

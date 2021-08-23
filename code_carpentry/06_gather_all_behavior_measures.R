weightx<-weight %>% 
  mutate(pair_status=status) %>% 
  select(-status) %>% 
  mutate(cohort=LETTERS[as.numeric(cohort)-104]) %>% 
  mutate(subjectID=paste(cohort,mouseID,sep="."))


dsrank<-dss.df.all %>% 
  mutate(subjectID=row.names(.),
         ds_rank=rank,
         cohort=as.character(cohort)) %>% 
  select(-rank) 


isirank<-m.isi %>% 
  as.data.frame %>% 
  mutate(isi_rank=row.names(.)) %>% 
  gather(cohort,mouseID,1:12) %>%  
  mutate(subjectID=paste(cohort,mouseID,sep="."))


glickorank<-lapply(df.glickos,function(x) x$ratings %>%
                     mutate(mouseID=Player,
                            glicko_rank=row.names(.)) %>% 
                     select(mouseID,Rating,glicko_rank,Win,Loss,Lag)) %>% 
  map2_df (.,names(.),~mutate(.x,cohort=.y)) %>% 
  mutate(subjectID=paste(cohort,mouseID,sep="."),
         mouseID=as.character(mouseID))


all<-weightx %>% 
  left_join(.,dsrank) %>% 
  left_join(.,isirank) %>% 
  left_join(.,glickorank) %>% 
  left_join(ind_cert_df) %>% 
  left_join(indiv_network_df)


all_behavior <- all %>% 
  left_join(.,df.sub_behav_perc %>% 
              mutate(mouseID = as.character(mouseID))) %>% 
  left_join(.,df.sub_behav_freq %>% 
              mutate(mouseID = as.character(mouseID)))


all_behavior %>% 
  group_by(cohort) %>% 
  mutate(win_prop = Win/sum(Win),
         loss_prop = Loss/sum(Loss)) %>% 
  ungroup() %>% 
  mutate(domgroup=ifelse(glicko_rank==1,"Alpha",ifelse(ds>=4.5,"Subdominant","Subordinate"))) -> all_behavior

all_behavior <- all_behavior %>%
  left_join(cluster)


all_behavior %>% 
  left_join(.,resultsdf %>% 
              select(cohort,gini.win,gini.lose,despotism), by = "cohort") -> all_behavior

colnames(all_behavior)

saveRDS(all_behavior,"results_statRDS/all_behavior.RDS")


#for multinomial analysis ======================================================
response_each<-df.groupx %>% map(get_response_freq) %>% 
  cohort_rbindlist() %>% 
  mutate(aggressorID=paste(cohort,aggressorID,sep=".")) %>% 
  select(- mouseID)


df1<-all_behavior %>% 
  mutate(domgroup=ifelse(glicko_rank==1,"Alpha",ifelse(ds>=4.5,"Subdominant","Subordinate"))) %>% 
  select(subjectID,glicko_rank,domgroup,ds)

df2<-all_behavior %>% 
  mutate(domgroup=ifelse(glicko_rank==1,"Alpha",ifelse(ds>=4.5,"Subdominant","Subordinate"))) %>% 
  select(subjectID,glicko_rank,domgroup,ds)
colnames(df2)<-c("aggressorID","agg_glicko_rank","agg_domgroup","ds")

response_each_rank<-response_each %>% 
  left_join(.,df1) %>% 
  left_join(.,df2) %>% 
  mutate(agg_glicko_rank=as.numeric(agg_glicko_rank),
         glicko_rank=as.numeric(glicko_rank),
         GD=as.numeric(GD),
         response=factor(response))


p1<-panel1 %>% 
  filter(is.na(diff_cd45)) %>% 
  mutate(cohort=ifelse(housing=="group",str_sub(sample_ID,2,4),NA)) %>% 
  mutate(ID= gsub('.*_(.*)','\\1',str_sub(sample_ID,1,-5))) %>% 
  mutate(ID= gsub('.*-(.*)','\\1',ID)) %>% 
  as.data.frame() 

p1x<-p1 %>% 
  mutate_at(.funs = funs(perc = round(./`Immune cells_count`*100,2)),
            .vars = vars(DCs_count:Neutro_count)) %>% 
  select_at(.vars= vars(batch,cohort,ID,housing,DCs_count_perc:Neutro_count_perc))


p2<-panel2 %>% 
  filter(is.na(diff_cd45)) %>% 
  mutate(cohort=ifelse(housing=="group",str_sub(sample_ID,2,4),NA)) %>% 
  mutate(ID= gsub('.*_(.*)','\\1',str_sub(sample_ID,1,-5))) %>% 
  mutate(ID= gsub('.*-(.*)','\\1',ID)) %>% 
  as.data.frame() 

p2x<-p2 %>% 
  mutate_at(.funs = funs(perc = round(./`Immune cells_count`*100,2)),
            .vars = vars(`B cells_count`:`Helper T_count`)) %>% 
  select_at(.vars= vars(batch,cohort,ID,housing, `B cells_count_perc`:`Helper T_count_perc`))


flow_df<-left_join(p1x,p2x) %>% 
  mutate(cohort=ifelse(housing=="group",LETTERS[as.numeric(cohort)-104],NA))

df_pair<-flow_df %>% 
  filter(housing=="pair") %>%
  mutate_at(.funs = funs(PD09 = round(.,2)),
            .vars = vars(`DCs_count_perc`:`Helper T_count_perc`)) %>% 
  mutate(pairID=str_sub(ID,-1),
         mouseID=str_sub(ID,1,-2)) %>% 
  select_at(.vars = vars(batch,pairID,mouseID,`DCs_count_perc_PD09`:`Helper T_count_perc_PD09`))


df_group<-flow_df %>% 
  filter(housing=="group") %>%
  mutate_at(.funs = funs(GD14 = round(.,2)),
            .vars = vars(`DCs_count_perc`:`Helper T_count_perc`)) %>% 
  mutate(mouseID=ID) %>% 
  select_at(.vars = vars(cohort,mouseID,`DCs_count_perc_GD14`:`Helper T_count_perc_GD14`))

  
flow<-all_behavior %>% 
  left_join(.,df_pair) %>%
  left_join(.,df_group)




## tidier flow data please
flow_dfx<-flow_df %>% 
  mutate(NLR=Neutro_count_perc/(`B cells_count_perc`+`T cells_count_perc`+`NK cells_count_perc`),
         MLR=Mono_count_perc/(`B cells_count_perc`+`T cells_count_perc`+`NK cells_count_perc`),
         CD4_CD8_ratio=`Helper T_count_perc`/`Cytotoxic T_count_perc`) 

flow_pair<-flow_dfx %>% 
  filter(housing=="pair") %>%
  mutate(cohort=ifelse(batch=="C","Pair-control",cohort)) %>% 
  gather(cell_type,perc,5:18) %>% 
  mutate(pairID=str_sub(ID,-1),
         mouseID=substr(ID,1,str_length(ID)-1)) %>% 
  select(-cohort,-ID) 

flow_pair<-flow_pair %>% left_join(.,all_behavior %>% 
                          select(batch,cohort,pairID,mouseID,pair_status,subjectID))




flow_group<-flow_dfx %>% 
  filter(housing=="group") %>% 
  gather(cell_type,perc,5:18) %>% 
  mutate(subjectID=paste(cohort,ID,sep=".")) %>% 
  select(-ID) %>% 
  left_join(.,all_behavior %>% 
              select(batch,cohort,pairID,mouseID,pair_status,subjectID))

flow_all<-rbind(flow_group,flow_pair) %>% 
  mutate(cohort=ifelse(is.na(cohort),"Pair-comparison",cohort)) %>% 
  mutate(subjectID=paste(cohort,mouseID,sep=".")) 


#I know this is not the tidiest code... :p

saveRDS(flow_all,'results_statRDS/flow_all.RDS')

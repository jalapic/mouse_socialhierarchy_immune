

all_behaviorx <- all_behavior %>% 
  left_join(cort_final) %>% 
  select(-time_out)

flowx<-flow_all %>% 
  filter(cohort !="Pair-comparison") %>% 
  mutate(timepoint = ifelse(housing == "group","GD14","PD09")) %>% 
  mutate(value=perc) %>% 
  mutate(measure = ifelse(cell_type == "CD4_CD8_ratio"|cell_type == "MLR"|cell_type =="NLR",
                            cell_type,
                            str_sub(cell_type,1,-12))) %>% 
  select(subjectID, cohort,timepoint,measure,value) 


fkbp5x<-fkbp5 %>% 
  filter(timepoint == "GD14") %>% 
  mutate(cohort = LETTERS[cohort - 104]) %>% 
  mutate(subjectID = as.character(glue('{cohort}.{mouseID}'))) %>% 
  mutate(CpG_mean = Mean) %>% 
  select(subjectID, cohort, timepoint, CpG_1,CpG_2,CpG_3,Mean) %>% 
  gather(measure,value,4:7)

measures <- rbind(flowx,fkbp5x)


measuresx <- measures %>% 
  mutate(measure = paste(measure,timepoint,sep = "_")) %>% 
  select(-timepoint) %>% 
  spread(measure,value)


alldata <- all_behaviorx %>% 
  left_join(.,measuresx) %>%  
  mutate(subjectID=paste(cohort,mouseID,sep="."))

colnames(alldata)


saveRDS(alldata,"results_statRDS/alldata.RDS")



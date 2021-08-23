colnames(rnaseq_rawcounts)[1] <- 'ensgene'
colnames(rnaseq_rawcounts) <- gsub('.trim.sam.counts','',colnames(rnaseq_rawcounts))

rnaseq_rawcounts %>% 
  select_if(!grepl('Spleen',names(.))) -> liver_counts


rnaseq_rawcounts %>% 
  select_if(!grepl('Liver',names(.))) -> spleen_counts


rnaseq_sampleid %>% 
  mutate(cohortx = as.numeric(str_sub(tissueID,1,3)),
         mouseID = as.numeric(str_sub(tissueID,5))) %>% 
  mutate(cohort = LETTERS[cohortx-104]) %>% 
  mutate(subjectID = glue('{cohort}.{mouseID}')) %>% 
  mutate(tissue = ifelse(grepl('Liver',sampleID), 'Liver','Spleen')) %>% 
  arrange(sampleID) %>% 
  select(sampleID, tissue, subjectID)-> rnaseq_sampleid 



rnaseq_sampleid %>% 
  filter(grepl('Liver',sampleID)) %>% 
  .$subjectID -> newcolnames_liver

rnaseq_sampleid %>% 
  filter(grepl('Spleen',sampleID)) %>% 
  .$subjectID -> newcolnames_spleen #yeah they are identical but I prefer this way

colnames(liver_counts)[2:length(liver_counts)] <- newcolnames_liver
colnames(spleen_counts)[2:length(spleen_counts)] <- newcolnames_spleen

colnames(liver_counts)
colnames(spleen_counts)



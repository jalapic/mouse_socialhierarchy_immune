## Behavior data
cohort105 <- read_csv("data_raw/behavior/csv/cohort105_wl_startend.csv")
cohort106 <- read_csv("data_raw/behavior/csv/cohort106_wl_startend.csv")
cohort107 <- read_csv("data_raw/behavior/csv/cohort107_wl_startend.csv")
cohort108 <- read_csv("data_raw/behavior/csv/cohort108_wl_startend.csv")
cohort109 <- read_csv("data_raw/behavior/csv/cohort109_wl_startend.csv")
cohort110 <- read_csv("data_raw/behavior/csv/cohort110_wl_startend.csv")
cohort111 <- read_csv("data_raw/behavior/csv/cohort111_wl_startend.csv")
cohort112 <- read_csv("data_raw/behavior/csv/cohort112_wl_startend.csv")
cohort113 <- read_csv("data_raw/behavior/csv/cohort113_wl_startend.csv")
cohort114 <- read_csv("data_raw/behavior/csv/cohort114_wl_startend.csv")
cohort115 <- read_csv("data_raw/behavior/csv/cohort115_wl_startend.csv")
cohort116 <- read_csv("data_raw/behavior/csv/cohort116_wl_startend.csv")



df.pair<- read_csv("data_raw/behavior/pair_status.csv") %>% 
  mutate(mouseID=as.character(mouseID))


# Bodyweight data
startbw <- read_csv("data_raw/bodyweight/start_bodyweight.csv") %>% 
  mutate(cohort=as.character(cohort),
         mouseID=as.character(mouseID))
endbw <- read_csv("data_raw/bodyweight/end_bodyweight.csv") 
# PD09bw<- read_csv("data_raw/bodyweight/PD09_bodyweight_blood_collection.csv")


weight<-left_join(startbw,endbw) %>% left_join(.,df.pair)

# write.csv(weight,"data_clean/bodyweight_clean.csv",row.names = F)
# Immunophenotyping data 
panel1 <- read_csv("data_raw/flow_cytometry_data/panel1_count_data_all.csv") 
panel2 <- read_csv("data_raw/flow_cytometry_data/panel2_count_data_all.csv") 


# fkbp5 methylation data 
fkbp5 <- read_csv("data_raw/pyrosequencing_data/fkbp5.csv") 


# CORT data

cort_names <- list.files(path ="data_raw/cort_data",pattern = '*_reads.csv')
cort_list <- list()

cort_list <- lapply(glue('data_raw/cort_data/{cort_names}'), read_csv)
names(cort_list) <- sub('\\_.*','',cort_names)

assay_id <- read_csv('data_raw/cort_data/assay_id_info.csv')

rnaseq_rawcounts <- read_csv('data_raw/rnaseq_data/counts_liver_spleen.csv') 


rnaseq_sampleid <- read_csv('data_raw/rnaseq_data/sample_id.csv') 

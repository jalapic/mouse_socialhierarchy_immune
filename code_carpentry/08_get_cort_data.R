# # install.packages('drc') # for fitting four-parameter logistic regression
# 
# 
# 
# #Dilution factor 1:125
# dilution_factor = 125
# 
# 
# # function to get rid of NA in Well ID column
# fill_well_id <- function(df){
#   for(i in 1:nrow(df)){
#     ifelse(is.na(df$`Well ID`[i]),
#            df$`Well ID`[i]<-df$`Well ID`[i-1],
#            df$`Well ID`[i])
#   }
#   return(df)
# }
# 
# # cleaning up
# cort_list <- lapply(cort_list, function(x) x %>%
#                       fill_well_id %>%
#                       select(-`Conc/Dil`) %>%
#                       rename(well_ID =`Well ID`,
#                              well_location = Well,
#                              OD = `450`))
# 
# lapply(cort_list,head)
# 
# 
# 
# # calculate %B/B0 for all reading first
# get_B_B0_perc <- function(df){
#   df %>%
#     filter(grepl('NSB', well_ID)) %>%
#     summarize(B0 = mean(OD)) %>%
#     unlist()-> NSB
# 
#   df1 <- df %>%
#     mutate(net_OD = OD - NSB)
# 
#   df1 %>%
#     filter(grepl('BLK', well_ID)) %>%
#     summarize(B0 = mean(net_OD)) %>%
#     unlist()-> B0
# 
#   df2 <- df1 %>%
#     mutate(B_B0_perc = net_OD/B0*100) %>%
#     filter(!grepl('NSB', well_ID)) %>%
#     filter(!grepl('BLK', well_ID))
# 
#   return(df2)
# }
# 
# cort_B_B0_perc_list <- cort_list %>% map(get_B_B0_perc)
# 
# 
# 
# # data frame listing concentration & well ID; same throughout batches
# conc_df = data.frame(
# 
#   Standard  = cort_list[[1]] %>%
#     filter(grepl('STD',well_ID)) %>%
#     select(well_ID) %>%
#     unique(),
# 
#   conc = 5000/(2^(0:8))
# )
# 
# conc_df
# 
# # Standard Curve for each data frame (each plate)
# 
# # According to the manual I have to fit 4PL (4-parameter logistic function)
# # stole code from: https://stats.stackexchange.com/questions/61144/how-to-do-4-parametric-regression-for-elisa-data-in-r
# 
# 
# fit_4PL_standard_curve <- function(df){
#   library('drc')
#     std_df <- df %>%
#     filter(grepl('STD',well_ID)) %>%
#     left_join(conc_df)
# 
#   std_mod <- drm(B_B0_perc ~ conc, data = std_df ,
#                  fct = LL.4(names=c("Slope","Lower Limit","Upper Limit","ED50")))
# 
#   return(std_mod)
# }
# 
# 
# std_model_list <- cort_B_B0_perc_list %>% map(fit_4PL_standard_curve)
# 
# # dev.off()
# for (i in 1:length(std_model_list)){
#   plot(std_model_list[[i]], main = names(std_model_list[i]))
# }
# 
# 
# for (i in 1:length(std_model_list)){
#   print(names(std_model_list)[i])
#   print(summary(std_model_list[[i]]))
# 
# }
# 
# 
# # once done inspecting model and confirmed each model is correct,
# 
# 
# get_X_4PL_curve <- function(model, Y) {
#   # Y = (A-D)/(1+(X/C)^B) + D, where Y = %B/B0, X = concentration
# 
#   A = model$coefficients[3] # A = Upper limit
#   D = model$coefficients[2] # D = Lower Limit
#   B = model$coefficients[1] # B = Slope
#   C = model$coefficients[4] # C = ED50
# 
#   # So, to get X based on Y,
# 
#   X = C*((A-D)/(Y-D)-1)^(1/B)
# 
#   return(X)
# }
# 
# 
# # I 'could have' written code wihtout forloop but oh well, :P
# cort_conc_list <- list()
# 
# for(i in 1:length(cort_B_B0_perc_list)){
# 
#   cort_conc_list[[i]] <- cort_B_B0_perc_list[[i]] %>%
#     filter(!grepl('STD',well_ID)) %>%
#     mutate(conc = get_X_4PL_curve(std_model_list[[i]],B_B0_perc))
#   names(cort_conc_list)[i] = names(cort_B_B0_perc_list[i])
# }
# 
# cort_conc_list %>% map(head) # sanity check
# 
# fill_name_id <- function(df){
#   for(i in 1:nrow(df)){
#     ifelse(is.na(df$`Name`[i]),
#            df$`Name`[i]<-df$`Name`[i-1],
#            df$`Name`[i])
#   }
#   return(df)
# }
# 
# cort_conc_list <- cort_conc_list %>% map(fill_name_id)
# 
# conc_df <- cort_conc_list %>% map2_dfr(.,names(.), ~mutate(.x, batchID = .y))
# 
# # don't forget to multiply by dilution factor!
# conc_df %>%
#   mutate(final_conc_ng_ul = conc * dilution_factor/1000) %>%
#   group_by(well_ID,  batchID) %>%
#   summarize(cort_ng_ul = mean(final_conc_ng_ul,na.rm = T)) %>%
#   select(well_ID, cort_ng_ul, batchID)  -> final_cort_conc
# 
# 
# 
# final_cort_conc %>%
#   ggplot(aes(x = cort_ng_ul))+
#   geom_histogram()+
#   facet_wrap(~batchID)+
#   theme_bw()
# 
# 
# final_cort_conc %>%
#   ggplot(aes(batchID, cort_ng_ul))+
#   geom_boxplot(alpha = 0.2, outlier.color = NA)+
#   geom_jitter(width = 0.1)+
#   theme_bw()
# 
# str_sub(conc_df$well_ID,4)
# # now match the ID and batch/well ID
# final_cort_conc %>%
#   mutate(assay_id = glue('{str_sub(batchID,6)}{str_sub(well_ID,4)}')) %>%
#   left_join(assay_id) %>%
#   ungroup() %>%
#   mutate(cort_timepoint = sub('\\_.*','',sample_info)) %>%
#   mutate(mouseID = sub('\\..*','',sub('.*\\_','',sample_info))) %>%
#   mutate(pairID = sub('.*\\.','',sample_info)) %>%
#   mutate(batch = sub('Batch','',batchID)) %>%
#   select(-well_ID,-assay_id, -batchID)-> temp1
# 
# head(temp1)
# 
# temp_pre <- temp1 %>%
#   filter(cort_timepoint == "Pre")
# 
# temp_post <- temp1 %>%
#   filter(cort_timepoint != "Pre")
# 
# 
# temp_pre
# temp_pre %>%
#   left_join(all_behavior %>%
#               select(mouseID, pairID, batch, cohort)) -> cort_pre
# 
# table(temp1$batch,temp1$mouseID)
# 
# table(temp_post$batch,temp_post$mouseID)
# 
# temp_post %>%
#   select(-batch) %>%
#   mutate(cohort = as.numeric(mouseID)) %>%
#   mutate(mouseID = pairID) %>%
#   mutate(cohort = LETTERS[cohort-104]) %>%
#   select(-pairID) %>%
#   left_join(all_behavior %>%
#               select(mouseID, cohort,batch,pairID)) -> cort_post
# 
# table(cort_post$cohort)
# 
# cort_final <- rbind(cort_pre,cort_post) %>%
#   select(-sample_info) %>%
#   arrange(batch,cohort) %>%
#   spread(cort_timepoint,cort_ng_ul) %>%
#   rename(cort_pre = Pre,
#          cort_post = Post)
# saveRDS(cort_final,"data_clean/cort_final.RDS")
# 
cort_final <- readRDS("data_clean/cort_final.RDS")
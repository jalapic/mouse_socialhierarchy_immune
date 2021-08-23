
# ME data frame
liver_MEs <- readRDS("results_RNAseqRDS/Liver_MEs.RDS") %>% 
  rownames_to_column("sampleID") %>% 
  left_join(rnaseq_sampleid) %>% 
  select(-tissue,-sampleID)
spleen_MEs <- readRDS('results_RNAseqRDS/Spleen_MEs.RDS') %>% 
  rownames_to_column("sampleID") %>% 
  left_join(rnaseq_sampleid) %>% 
  select(-tissue,-sampleID)

colnames(liver_MEs) <- gsub('ME','Liver_',colnames(liver_MEs))
colnames(spleen_MEs) <- gsub('ME','Spleen_',colnames(spleen_MEs))

liver_MEs %>% left_join(spleen_MEs) -> all_MEs


# trait data frame 
alldata_MEs <- all_MEs %>% 
  left_join(alldata %>%  
    dplyr::select(subjectID, 
        GD01_BW, GD14_bw, spleen_weight_mg, AGD, pair_status, domgroup, 
        cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
        `B cells_GD14`, `CD4_CD8_ratio_GD14`, `Cytotoxic T_GD14`, `DCs_GD14`,
        `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
        `Neutro_GD14`, `NK cells_GD14`, `NLR_GD14`, `T cells_GD14`,
        `Lymphoid DCs_GD14`,`Myeloid_DCs_GD14`,
        `B cells_PD09`, `CD4_CD8_ratio_PD09`, `Cytotoxic T_PD09`, `DCs_PD09`,
        `Helper T_PD09`, `Macro_PD09`, `MLR_PD09`, `Mono_PD09`,
        `Neutro_PD09`, `NK cells_PD09`, `NLR_PD09`, `T cells_PD09`,
        `Lymphoid DCs_PD09`,`Myeloid_DCs_PD09`
  ) %>% 
    mutate_if(is.numeric,scale) %>%
    mutate(status = ifelse(domgroup == "Alpha", 1,-1)) %>% 
    select(-domgroup) %>% 
    mutate(pair_status = ifelse(pair_status == "dom", 1,-1)))

# library("Hmisc")
# res2 <- rcorr(as.matrix(all_MEs %>% select(-subjectID)))
# corrplot::corrplot(res2$r, type="upper",
#          p.mat = res2$P, sig.level = 0.01, insig = "blank") # correlation is bs

colnames(alldata_MEs)

alldata_MEs %>% 
  gather(key, value, 27:65) -> ME_trait_df

ME_trait_df %>% 
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~ key, scale = "free")

ME_trait_df %>% 
  ggplot(aes(value,Liver_pink))+
  geom_point(shape = 21)+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~ key, scale = "free")

ME_trait_df %>% 
  ggplot(aes(value,Liver_blue))+
  geom_point(shape = 21)+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~ key, scale = "free")


ME_trait_df %>% 
  split(.$key) -> ME_trait_list


ME_trait_list[33] -> x
str(x)
x %>% as.data.frame


# linear regression
do_lm <- function(x){
x %>% bind_rows() -> x
  x$key %>% unique() -> my_name

x %>% 
  select(-subjectID, -key) -> xx


summary(lm(xx[,1] ~ xx$value))$coefficients[2,] %>% 
  bind_rows() %>% 
  mutate(trait = my_name) -> lm_result

  for(i in 2:(length(xx)-1)){
    summary(lm(xx[,i] ~ xx$value))$coefficients[2,] %>% 
      bind_rows() %>% 
      mutate(trait = my_name) -> temp
    rbind(lm_result, temp) -> lm_result
  }
lm_result$module_name <- colnames(xx)[1:(length(xx)-1)]
  return(lm_result)
}


lapply(ME_trait_list,do_lm) -> lm_result_list
lm_result_list$AGD %>% head(43) %>% as.data.frame()

saveRDS(lm_result_list, "results_RNAseqRDS/lm_result_list.RDS")


# spearman rank correlation
do_spearman_corr <- function(x){
  x %>% bind_rows() -> x
  x$key %>% unique() -> my_name
  
  x %>% 
    select(-subjectID, -key) -> xx
  
  temp <- cor.test(xx[,1],xx$value,method ="spearman", use = "complete.obs",exact=FALSE)
  str(temp)
  temp$p.value -> p.val
  temp$estimate -> rho
  cbind(p.val, rho) %>% 
    as.data.frame() %>% 
  mutate(trait = my_name) -> lm_result
  
  for(i in 2:(length(xx)-1)){
    temp <- cor.test(xx[,i],xx$value,method ="spearman", use = "complete.obs",exact=FALSE)
    str(temp)
    temp$p.value -> p.val
    temp$estimate -> rho
    cbind(p.val, rho) %>% 
      as.data.frame() %>% 
      mutate(trait = my_name) -> temp
    rbind(lm_result, temp) -> lm_result
  }
  lm_result$module_name <- colnames(xx)[1:(length(xx)-1)]
  return(lm_result)
}


lapply(ME_trait_list,do_spearman_corr) -> sp_corr_result_list

saveRDS(sp_corr_result_list, "results_RNAseqRDS/sp_corr_result_list.RDS")




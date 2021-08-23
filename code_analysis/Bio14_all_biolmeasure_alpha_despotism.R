

alphadata <- alldata %>% 
  filter(glicko_rank == 1)

alphadata %>% 
 select(subjectID,despotism,win_prop,
       "GD01_BW", "GD14_bw", "spleen_weight_mg", "AGD",
        "B cells_GD14", "B cells_PD09", "CD4_CD8_ratio_GD14", "CD4_CD8_ratio_PD09",
        "CpG_1_GD14", "CpG_2_GD14", "CpG_3_GD14", 
        "Cytotoxic T_GD14", "Cytotoxic T_PD09", "Helper T_GD14", "Helper T_PD09",
        "Lymphoid DCs_GD14", "Lymphoid DCs_PD09", "Macro_GD14", "Macro_PD09",
        "MLR_GD14", "MLR_PD09" , "Mono_GD14", "Mono_PD09",
        "Myeloid_DCs_GD14", "Myeloid_DCs_PD09", "Neutro_GD14", "Neutro_PD09",
        "NK cells_GD14", "NK cells_PD09", "T cells_GD14", "T cells_PD09") %>% 
  gather(bio,value,4:34) -> alphadf


alphadf %>% 
  ggplot(aes(despotism,value))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap( ~bio, scales = "free_y")



# STAT 
alphadf %>% 
  split(.$bio) -> alpha_list
result_list_despotism <- list()

for (i in 1: length(alpha_list))
{
  result_list_despotism[[i]]<-
    brm(value ~ despotism,
        data = alpha_list[[i]],
        file = glue('results/despotism_{names(alpha_list)[[i]]}.Rds')
    )
}

names(result_list_despotism) <- names(alpha_list)

lapply(result_list_despotism, 
       function(x) x %>% 
         spread_draws(b_despotism) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_despotism



plotdf_despotism %>% 
  filter(.width == 0.95) %>% 
  as.data.frame()

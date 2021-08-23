
# spleen weight and post cort? 

colnames(alldata)

alldata %>% 
  filter(cort_post>5) %>% 
  gather(fkbp5_site,fkbp5_value, c(CpG_1_GD14,CpG_2_GD14, CpG_3_GD14)) %>% 
  filter(fkbp5_value >50) %>% 
  unique() -> df

table(df$cohort)
table(alldata$cohort)

alldata %>% 
  ggplot(aes(spleen_weight_mg, cort_pre, color = domgroup))+
  geom_point(size = 4, alpha = 0.5)

alldata %>% 
  ggplot(aes(spleen_weight_mg, cort_post, color = domgroup))+
  geom_point(size = 4, alpha = 0.5)

df %>% 
  ggplot(aes(fkbp5_value, cort_pre, color = spleen_weight_mg))+
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~fkbp5_site, scale = "free")+
  scale_color_viridis()



df %>% 
  ggplot(aes(fkbp5_value, cort_post, color = spleen_weight_mg))+
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~fkbp5_site, scale = "free")+
  scale_color_viridis()




df %>% 
  ggplot(aes(fkbp5_value, cort_post-cort_pre, color = spleen_weight_mg))+
  geom_point(size = 4, alpha = 0.5) +
  facet_wrap(~fkbp5_site, scale = "free")+
  # geom_smooth(method = "lm")+
  scale_color_viridis() # meh I am not convinced 


df %>%
  ggplot(aes(fkbp5_value))+
  geom_histogram()+
  facet_wrap(~fkbp5_site)





#cORRELATION HEATMAP (ALL, BY DOMGGROUP - so total four columns  )


cormat_df <- 
  alldata %>% 
  mutate(status = ifelse(domgroup == "Alpha", 1,0)) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", 1,0)) %>% 
  column_to_rownames('subjectID') %>% 
  select(
    GD01_BW, GD14_bw, spleen_weight_mg, AGD, pair_status,
    status,ds, glicko_rank, win_prop, loss_prop, ind_cert,
    preem_flee_perc, flee_perc, sub_pos_perc, freeze_perc, 
    despotism, cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
    `B cells_GD14`, CD4_CD8_ratio_GD14, `Cytotoxic T_GD14`, `DCs_GD14`,
    `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
    `Neutro_GD14`, `NK cells_GD14`, `T cells_GD14`
  )

cormat_df


library("Hmisc")
res2 <- rcorr(as.matrix(cormat_df))
corrplot::corrplot(res2$r, type="upper", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank") # correlation is bs


library("PerformanceAnalytics")
chart.Correlation(cormat_df %>% 
                    select(ds, cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
                           `B cells_GD14`, CD4_CD8_ratio_GD14, `Cytotoxic T_GD14`, `DCs_GD14`,
                           `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
                           `Neutro_GD14`, `NK cells_GD14`, `T cells_GD14`),
                  histogram=TRUE, pch=19) # welp, not that helpful 



# only subordiantes ==================================================

alldata %>% 
  filter(domgroup == "Subordinate") -> sub

sub %>% 
  ggplot(aes(loss_prop))+
  geom_histogram()


# win loss proportion

png(filename = glue("results_figures/CORT_lossprop_onlysub.png"),
    width = 6, height = 5, units = "cm", res = 600)
sub %>% 
  ggplot(aes(loss_prop,cort_post))+
  geom_point(shape = 21)+
  theme_bw(base_size = 6)+
  labs(x = "Loss proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only")

invisible(dev.off())


png(filename = glue("results_figures/CORT_winprop_onlysub.png"),
    width = 6, height = 5, units = "cm", res = 600)
sub %>% 
  ggplot(aes(win_prop,cort_post))+
  geom_point(shape = 21)+
  theme_bw(base_size = 6)+
  labs(x = "Win proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only")

invisible(dev.off())


# Actual win and loss numbers 

png(filename = glue("results_figures/CORT_loss_onlysub.png"),
    width = 6, height = 5, units = "cm", res = 600)
sub %>% 
  ggplot(aes(Loss,cort_post))+
  geom_point(shape = 21)+
  theme_bw(base_size = 6)+
  labs(x = "Total loss observed", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only")

invisible(dev.off())


png(filename = glue("results_figures/CORT_win_onlysub.png"),
    width = 6, height = 5, units = "cm", res = 600)
sub %>% 
  ggplot(aes(Win,cort_post))+
  geom_point(shape = 21)+
  theme_bw(base_size = 6)+
  labs(x = "Total win observed", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only")

invisible(dev.off())


# Want to add pair housing social status 

png(filename = glue("results_figures/CORT_lossprop_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(loss_prop,cort_post, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  scale_fill_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Loss proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())


png(filename = glue("results_figures/CORT_wibprop_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(win_prop,cort_post, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  scale_fill_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Win proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())

# what about monocyte (you can see I am desperate)

png(filename = glue("results_figures/CORT_monocyte_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(Mono_GD14,cort_post, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Monocyte proportion (GD14)", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())


png(filename = glue("results_figures/CORT_monopd09_onlysub_pairstatus.png"),
    width = 6, height = 6, units = "cm", res = 600)
sub %>% 
  ggplot(aes(Mono_PD09,cort_post, color = pair_status))+
  geom_point(shape = 21)+
  scale_color_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  scale_fill_manual(values = c("purple4","orange"),label = c("Dominant","Subordinate"))+
  theme_bw(base_size = 6)+
  labs(x = "Monocyte proportion (PD09)",
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only",
       color = "Pair housing status",
       fill = "Pair housing status")+
  theme(legend.position = "top")

invisible(dev.off())




# ahhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

png(filename = glue("results_figures/CORT_winprop_onlysub_smooth.png"),
    width = 6, height = 5, units = "cm", res = 600)
sub %>% 
  ggplot(aes(win_prop,cort_post, fill = ds, color = ds))+
  geom_point(shape = 21)+
  geom_smooth(method = "lm", se = F)+
  theme_bw(base_size = 6)+
  scale_color_distiller(palette = "Spectral")+
  scale_fill_distiller(palette = "Spectral")+
  labs(x = "Win proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       title = "Subordinate mice only")

invisible(dev.off())


  png(filename = glue("results_figures/CORT_winprop.png"),
    width = 8, height = 5, units = "cm", res = 600)
df %>% 
  ggplot(aes(win_prop,cort_post, shape = km3_cluster, color = km3_cluster))+
  geom_point(fill = NA, alpha = 0.2)+
  geom_smooth(method = "lm")+
  theme_bw(base_size = 6)+
  labs(x = "Win proportion", 
       y = "Corticosterone (GD14) (ng/ml)",
       shape = "",
       color = "")+
  theme(legend.position = "top")

invisible(dev.off())





# what is going on with some high cort subordinates???


sub %>%
  filter(cort_post >400) %>%
  arrange(desc(cort_post)) %>% 
  select(subjectID, Lag, ind_cert, cort_post,glicko_rank, pair_status, km3_cluster)

















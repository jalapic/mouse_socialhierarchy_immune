my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}_second.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}_second.RDS"))}


limma_list$status %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_status = logFC,
         pval_status = P.Value) -> status_df

limma_list$cort %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_cort = logFC,
         pval_cort = P.Value) -> cort_df

df <- full_join(status_df, cort_df) %>% 
  mutate(Sig = ifelse(pval_status >= 0.05 & pval_cort >= 0.05, "None", 
                      ifelse(pval_status < 0.05 & pval_cort < 0.05,"Both",
                             ifelse(pval_status < 0.05,"Status-specific","CORT-specific")))) %>% 
  mutate(Sig = factor(Sig, levels =c("Both","Status-specific","CORT-specific","None")))

lim = 1.5

# png(filename = glue("results_figures/twoaxis_{my_tissue}_allgene.png"),
#     width = 7, height = 7.6, units = "cm", res = 600)

df %>% 
  ggplot(aes(logFC_status, logFC_cort))+
  geom_abline(slope = 1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_abline(slope = -1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_hline(yintercept = 0, color ='grey', linetype = 'dashed')+
  geom_vline(xintercept = 0, color ='grey', linetype = 'dashed')+
  geom_point(shape = 21, size = 0.3)+
  labs(x = "Status effect (Sub <-> Dom)",
       y = "CORT effect (Low <-> High)",
       color = "",
       title = glue("{my_tissue}"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "top")+
  scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))

# invisible(dev.off())



png(filename = glue("results_figures/twoaxis_{my_tissue}_allgene2z.png"),
    width = 7, height = 7.6, units = "cm", res = 600)


ggplot(df,aes(color = Sig))+
  geom_abline(slope = 1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_abline(slope = -1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_hline(yintercept = 0, color ='grey', linetype = 'dashed')+
  geom_vline(xintercept = 0, color ='grey', linetype = 'dashed')+
  geom_point(data = df %>% 
               filter(Sig == "None"), 
             aes(logFC_status, logFC_cort),
             color = "grey", #fill = "grey",
             shape = 21, size = 0.3)+
  geom_point(data = df %>% 
               filter(Sig == "Status-specific"), 
             aes(logFC_status, logFC_cort),
             color = "#E7B800",#fill = "#d8b365",
             shape = 21, size = 0.3)+
  geom_point(data = df %>% 
               filter(Sig == "CORT-specific"), 
             aes(logFC_status, logFC_cort),
             color = "#00AFBB" ,#fill = "#5ab4ac",
             shape = 21, size = 0.3)+
  geom_point(data = df %>% 
               filter(Sig == "Both"), 
             aes(logFC_status, logFC_cort),
             color = "#FC4E07",#fill  = "#fdc086",
             shape = 21, size = 0.3)+
  labs(x = "Status effect (Sub <-> Dom)",
       y = "CORT effect (Low <-> High)",
       color = "",
       title = glue("{my_tissue}"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "top")+
  scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))

invisible(dev.off())


# png(filename = glue("results_figures/twoaxis_{my_tissue}_allgene3.png"),
#     width = 7, height = 7.6, units = "cm", res = 600)

df %>% 
  ggplot(aes(logFC_status, logFC_cort, color = Sig))+
  geom_abline(slope = 1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_abline(slope = -1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_hline(yintercept = 0, color ='grey', linetype = 'dashed')+
  geom_vline(xintercept = 0, color ='grey', linetype = 'dashed')+
  geom_point(shape = 21, size = 3.3)+
  labs(x = "Status effect (Sub <-> Dom)",
       y = "CORT effect (Low <-> High)",
       color = "P-value < 0.05",
       fill = "P-value < 0.05",
       title = glue("{my_tissue}"))+
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.8,0.2),
        legend.key.height = unit(0,"mm"))+
  scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
  scale_color_manual(values = c("#FC4E07","#E7B800","#00AFBB" ,"grey"))-> ppp
# dev.off()


png(filename = glue("results_figures/legend_allgene.png"),
    width = 7, height = 7.6, units = "cm", res = 600)

leg <- get_legend(ppp)
grid.arrange(leg)

dev.off()


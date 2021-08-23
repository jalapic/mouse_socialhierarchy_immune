my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){
  decon_df <- read_csv("results_CIBERSORT/liver_immucc_expression_result.csv")  
}else{decon_df <- read_csv("results_CIBERSORT/Spleen_immucc_expression_result.csv")  }

colnames(decon_df)[1]<- "subjectID"
decon_df <- decon_df[,-c(10:12)]
# head(decon_df)
# decon_df %>% View()

png(filename = glue("results_figures/Deconvolution_{my_tissue}.png"),
    width = 18, height = 12, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  ggplot(aes(domgroup, value))+
  geom_line(aes(group = cohort), color = 'grey')+
  geom_point(aes(color = domgroup),alpha = 0.3, size = 2)+
  scale_color_manual(values = c("purple4", "orange"))+
    facet_wrap(~key, nrow = 2)+
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))+
  theme(legend.position = "none") 

invisible(dev.off())


png(filename = glue("results_figures/Deconvolution_{my_tissue}_boxjitter.png"),
    width = 18, height = 12, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  ggplot(aes(domgroup, value, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  facet_wrap(~key, nrow = 2)+
  scale_color_manual(values = c("purple4", "orange"))+
  scale_fill_manual(values = c("purple4", "orange"))+
  theme(legend.position = "none") +
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))

invisible(dev.off())

 # stacked bar ============================================================
read_csv("results_CIBERSORT/liver_immucc_expression_result.csv")[,1:9]  %>%
  rename(subjectID = Mixture) %>% 
  summarize_all(.funs = mean) %>% 
  mutate(subjectID = "Liver") %>% 
  gather(key,value, 2:ncol(.)) -> liver_df

read_csv("results_CIBERSORT/spleen_immucc_expression_result.csv")[,1:9]  %>%
  rename(subjectID = Mixture) %>% 
  summarize_all(.funs = mean) %>% 
  mutate(subjectID = "Spleen") %>% 
  gather(key,value, 2:ncol(.)) -> spleen_df



png(filename = glue("results_figures/Deconvolution_stacked.png"),
    width = 8, height = 10, units = "cm", res = 600)

rbind(liver_df, spleen_df) %>% 
  ggplot(aes(x = subjectID, y =value, fill = key))+
  # geom_point(stat='identity')+
  geom_col(position = "stack",color ='white', alpha = 0.8)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_minimal(base_size = 10)+
  labs(x = "",
       y = "",
       fill = "Immune cell type")+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(vjust = 5, size = 10))
invisible(dev.off())


png(filename = glue("results_figures/Deconvolution_stacked_flip.png"),
    width = 14, height = 6, units = "cm", res = 600)

rbind(liver_df, spleen_df) %>% 
  ggplot(aes(x = reorder(subjectID,desc(subjectID)), y =value, fill = key))+
  # geom_point(stat='identity')+
  geom_col(position = "stack",color ='white', alpha = 0.8)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_minimal(base_size = 10)+
  # coord_flip()+
  labs(x = "",
       y = "",
       fill = "Immune cell type")+
  theme(panel.grid = element_blank(),
        # legend.position = "top", # for flip
        axis.text.y = element_text(vjust = 5, size = 10)) #change to hjust = 1.2 for coordflip
invisible(dev.off())



# Liver =======================================================================
my_tissue = "Liver"


if(my_tissue == "Liver"){
  decon_df <- read_csv("results_CIBERSORT/liver_immucc_expression_result.csv")  
}else{decon_df <- read_csv("results_CIBERSORT/Spleen_immucc_expression_result.csv")  }

colnames(decon_df)[1]<- "subjectID"
decon_df <- decon_df[,-c(10:12)]
# head(decon_df)
# decon_df %>% View()

png(filename = glue("results_figures/Deconvolution_{my_tissue}2.png"),
    width = 18, height = 8, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  filter(key %in% c("Macrophage", "Monocyte","B cell", "T cell")) %>% 
  mutate(key = factor(key, levels = c("Macrophage", "Monocyte","B cell", "T cell"))) %>% 
  ggplot(aes(domgroup, value))+
  geom_line(aes(group = cohort), color = 'grey')+
  geom_point(aes(color = domgroup),alpha = 0.3, size = 2)+
  scale_color_manual(values = c("purple4", "orange"))+
  facet_wrap(~key, nrow = 1, scales = "free_y")+
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))+
  theme(legend.position = "none") 

invisible(dev.off())


png(filename = glue("results_figures/Deconvolution_{my_tissue}_boxjitter2.png"),
    width = 18, height = 8, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  filter(key %in% c("Macrophage", "Monocyte","B cell", "T cell")) %>% 
  mutate(key = factor(key, levels = c("Macrophage", "Monocyte","B cell", "T cell"))) %>% 
  ggplot(aes(domgroup, value, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.00, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  facet_wrap(~key, nrow = 1, scales = "free_y")+
  scale_color_manual(values = c("purple4", "orange"))+
  scale_fill_manual(values = c("purple4", "orange"))+
  theme(legend.position = "none") +
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))

invisible(dev.off())


# Spleen ============================================================================
my_tissue = "Spleen"

if(my_tissue == "Liver"){
  decon_df <- read_csv("results_CIBERSORT/liver_immucc_expression_result.csv")  
}else{decon_df <- read_csv("results_CIBERSORT/Spleen_immucc_expression_result.csv")  }

colnames(decon_df)[1]<- "subjectID"
decon_df <- decon_df[,-c(10:12)]
# head(decon_df)
# decon_df %>% View()

png(filename = glue("results_figures/Deconvolution_{my_tissue}2.png"),
    width = 18, height = 8, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  filter(key %in% c("Macrophage", "Monocyte", "Neutrophil","B cell")) %>% 
  mutate(key = factor(key, levels = c("Macrophage", "Monocyte", "Neutrophil","B cell"))) %>% 
  ggplot(aes(domgroup, value))+
  geom_line(aes(group = cohort), color = 'grey')+
  geom_point(aes(color = domgroup),alpha = 0.3, size = 2)+
  scale_color_manual(values = c("purple4", "orange"))+
  facet_wrap(~key, nrow = 1, scales = "free_y")+
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))+
  theme(legend.position = "none") 

invisible(dev.off())


png(filename = glue("results_figures/Deconvolution_{my_tissue}_boxjitter2.png"),
    width = 18, height = 8, units = "cm", res = 600)


decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  filter(key %in% c("Macrophage", "Monocyte", "Neutrophil","B cell")) %>% 
  mutate(key = factor(key, levels = c("Macrophage", "Monocyte", "Neutrophil","B cell"))) %>% 
  ggplot(aes(domgroup, value, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.00, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  facet_wrap(~key, nrow = 1, scales = "free_y")+
  scale_color_manual(values = c("purple4", "orange"))+
  scale_fill_manual(values = c("purple4", "orange"))+
  theme(legend.position = "none") +
  labs(x = "Social status",
       y = "Estimated cell proportion",
       title = glue("{my_tissue}"))

invisible(dev.off())


# Now stat and posterior plots ============================================================

my_tissue = "Liver"
my_key = c("Macrophage", "Monocyte","B cell", "T cell")

my_tissue = "Spleen"
my_key = c("Macrophage", "Monocyte", "Neutrophil","B cell")

if(my_tissue == "Liver"){
  decon_df <- read_csv("results_CIBERSORT/liver_immucc_expression_result.csv")  
}else{decon_df <- read_csv("results_CIBERSORT/Spleen_immucc_expression_result.csv")  }

colnames(decon_df)[1]<- "subjectID"
decon_df <- decon_df[,-c(10:12)]

decon_df %>% 
  left_join(alldata %>% select(subjectID, domgroup,cohort)) %>%
  gather(key, value, 2:9) %>% 
  filter(domgroup != 'Subdominant') %>% 
  filter(key %in% my_key) %>% 
  mutate(key = factor(key, levels = my_key)) %>% 
  split(.$key) -> decon_list


decon_brms_results <- vector('list',length = length(my_key))

for (i in 1:length(my_key)){
  
  decon_list[[i]] -> df
  gsub(" ", "", names(decon_list)[[i]]) -> keyname
  names(decon_brms_results)[[i]] <- names(decon_list)[[i]] 
  decon_brms_results[[i]] <-  brm(value ~ domgroup,
      data = df,
      file = glue("results_RNAseqRDS/brms_{my_tissue}_{keyname}.RDS"))
  
}

lapply(decon_brms_results, 
       function(x) x %>% 
         spread_draws(b_domgroupSubordinate) %>% #change
         median_qi(.width = c(.95, .66))) %>% 
  map2_df(.,names(.), ~ mutate(.x, cell_type = .y)) -> plotdf_decon

plotdf_decon[,1:3] <- (-1)*plotdf_decon[,1:3]
colnames(plotdf_decon)[2:3] <- c( ".upper", ".lower")

min(plotdf_decon$.lower)*3 -> min
# 2.5 for spleen, 1.8 for liver
plotdf_decon %>% 
  rename(myx = b_domgroupSubordinate) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_decon) %>% 
  rename(myx = b_domgroupSubordinate) -> plotdf1

plotdf1 %>% 
  mutate(cell_type = factor(cell_type, levels = rev(my_key))) %>% 
  ggplot(aes(y = cell_type, x = myx ,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf1 %>% filter(.width == 0.66),
               aes(y = cell_type, yend = cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf1 %>% filter(.width == 0.95),
               aes(y = cell_type, yend = cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base()+
  labs(x = "Subordinate       --       Estimates (median ± 95% CI)       --       Dominant",
       y = "")+
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title.x = element_text(size = 10),
        plot.background = element_blank())+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf1 %>% filter(sig == "Sig"),
            aes(label = statsummary, x= min), hjust = "inward", color = "red", fontface = "bold")+ # hash it for liver!
  geom_text(data = plotdf1 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= min), hjust = "inward", color = "black", )+ #
  ggtitle(glue("{my_tissue}: Immune cell proportion difference")) -> plot1

print(plot1)


png(filename = glue("results_figures/deconvolution_brms_{my_tissue}_narrow.png"),
    width = 15, height = 5, units = "cm", res = 600)
print(plot1)
invisible(dev.off())



theme_set(theme_bw(base_size = 7))
















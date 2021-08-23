

genecat = "Hallmark"
# msigdbr(species = "Mus musculus", category ="H") -> hallmark_sets
# saveRDS(hallmark_sets,"results_RNAseqRDS/hallmark_sets.RDS")
my_gs_sets <- readRDS("results_RNAseqRDS/hallmark_sets.RDS")



my_tissue = "Liver"
my_tissue = "Spleen"

limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
limma_vdl <- readRDS(glue("results_RNAseqRDS/limma_vdl_{my_tissue}.RDS"))


my_gs_sets$gs_name %>% unique() -> my_gs_names

i = 36

my_hallmark <- gsub(glue("{genecat}_"),"",my_gs_names[i],ignore.case = T)
hallmark <- my_gs_names[i]
my_gs_sets %>% 
  filter(gs_name == hallmark) %>% 
  .$gene_symbol -> selected_genes


limma_list$status %>% 
  filter(symbol %in% selected_genes) %>% 
  select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>0,"Alpha genes","Subordinate genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("Alpha genes","Subordinate genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> deg_df

# count number of DEG in hallmark gene set 
deg_df %>% 
  group_by(Sig) %>% 
  count()

alldata %>% filter(domgroup == "Alpha") %>% .$subjectID -> alpha_id
alldata %>% filter(domgroup != "Alpha") %>% .$subjectID -> sub_id

limma_vdl$E %>% 
  as.data.frame() %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  filter(symbol %in% selected_genes) %>% 
  mutate(Alpha_Mean = rowMeans(select_if(., colnames(.) %in% alpha_id))) %>% 
  mutate(Sub_Mean = rowMeans(select_if(., colnames(.) %in% sub_id))) %>% 
  mutate(Mean = rowMeans(select(., -ensgene, -symbol))) %>% 
  left_join(deg_df)  %>% 
  select(symbol,Alpha_Mean,Sub_Mean,Mean, logFC,Sig, P.Value) -> hallmark_all_df



# plot 1 - easy
hallmark_all_df %>% 
  filter(!is.na(Sig)) %>% 
  arrange(desc(Sig), desc(Mean)) %>% 
  mutate(neworder = row_number()) %>% 
  ggplot(aes(reorder(symbol,-neworder), Mean,color = Sig))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "top")+
  scale_color_manual(values = c("purple4","orange","grey"))+
  labs(x = "",
       y = "Mean expression level",
       color = "",
       title = glue("{my_tissue}: {my_hallmark}"))
  
# plot 2 - divide into alpha and sub 
png(filename = glue("results_figures/hallmark_expressiondiff_{my_tissue}_{my_hallmark}.png"),
    width = 12, height = 24, units = "cm", res = 600)
hallmark_all_df %>% 
  filter(!is.na(Sig)) %>% 
  arrange(desc(Sig), Mean) %>% 
  mutate(neworder = nrow(.)-row_number()) %>% 
  ggplot(aes(y=reorder(symbol,-neworder), yend=reorder(symbol,-neworder),
             x=Alpha_Mean,xend = Sub_Mean, 
             color = Sig))+
  geom_segment()+
  theme(legend.position = "top")+
  scale_color_manual(values = c("purple4","orange","grey"))+
  scale_x_continuous(sec.axis = dup_axis(),expand=expansion(mult=c(0,0.0)),
                     breaks = c(0,2,4,6,8,10,12,14))+
  labs(x = "Mean expression level",
       y = "",
       color = "",
       title = glue("{my_tissue}: {my_hallmark}"))

dev.off()


## now for all 50 gene sets, what is the mean epxression level? 
alldata %>% filter(domgroup == "Alpha") %>% .$subjectID -> alpha_id
alldata %>% filter(domgroup != "Alpha") %>% .$subjectID -> sub_id

hallmark_expression_list <- vector('list', length = length(my_gs_names))

for (i in 1: length(my_gs_names)){
  my_hallmark <- gsub(glue("{genecat}_"),"",my_gs_names[i],ignore.case = T)
  hallmark <- my_gs_names[i]
  my_gs_sets %>% 
    filter(gs_name == hallmark) %>% 
    .$gene_symbol -> selected_genes
  
  
  limma_list$status %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                        ifelse(logFC>0,"Alpha genes","Subordinate genes"))) %>%
    mutate(Sig = factor(Sig, levels = c("Alpha genes","Subordinate genes", "N.S."))) %>% 
    mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
    unique() -> deg_df
  
  limma_vdl$E %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    left_join(grcm38 %>% select(ensgene, symbol)) %>% 
    filter(symbol %in% selected_genes) %>% 
    mutate(Alpha_Mean = rowMeans(select_if(., colnames(.) %in% alpha_id))) %>% 
    mutate(Sub_Mean = rowMeans(select_if(., colnames(.) %in% sub_id))) %>% 
    mutate(Mean = rowMeans(select(., -ensgene, -symbol))) %>% 
    left_join(deg_df)  %>% 
    select(symbol,Alpha_Mean,Sub_Mean,Mean, logFC,Sig, P.Value) -> hallmark_all_df
  
  hallmark_expression_list[[i]] <- hallmark_all_df %>% 
    mutate(Hallmark = my_hallmark)
  
}

rbindlist(hallmark_expression_list) -> big_plotdf


png(filename = glue("results_figures/hallmark_expressionlevel_{my_tissue}.png"),
    width = 12, height = 24, units = "cm", res = 600)

big_plotdf %>% 
  filter(!is.na(Sig)) %>% 
  filter(Mean > 0) %>% 
  ggplot(aes(Mean,reorder(Hallmark,desc(Hallmark)),color = Sig, fill = Sig))+
  geom_jitter(shape = 21,alpha = 0.5,width = 0.1, height = 0.2,size=0.8)+
  theme_bw(base_size = 10)+
  theme(legend.position = "top")+
  scale_x_continuous(sec.axis = dup_axis(),expand=expansion(mult=c(0,0.0)),
                     breaks = c(0,2,4,6,8,10,12,14))+
  labs(x = "Mean expression level",
       y = "",
       color = "",
       fill = "",
       title = glue("{my_tissue}"))+
  scale_color_manual(values = c("purple4","orange","grey"))+
  scale_fill_manual(values = c("purple4","orange","grey"))


invisible(dev.off())


# now number of DEGs Alpha, Sub, N.S. for each hallmark gene set
# perhaps stacked bar - proportion
png(filename = glue("results_figures/hallmark_DEGproportion_{my_tissue}.png"),
    width = 12, height = 24, units = "cm", res = 600)

big_plotdf %>% 
  filter(!is.na(Sig)) %>% 
  # filter(Mean > 0) %>% 
  group_by(Hallmark,Sig) %>% 
  count() %>% 
  ggplot(aes(n,reorder(Hallmark,desc(Hallmark)),color = Sig, fill = Sig))+
  geom_col(position = "fill",color = NA)+
  theme_bw(base_size = 10)+
  theme(legend.position = "top")+
    labs(x = "Mean expression level",
       y = "",
       color = "",
       fill = "",
       title = glue("{my_tissue}"))+
  scale_x_continuous(sec.axis = dup_axis(),expand=expansion(mult=c(0,0.0)),
                     breaks = c(0,0.5,1))+
  scale_color_manual(values = c("purple4","orange","grey"))+
  scale_fill_manual(values = c("purple4","orange","grey"))
  
invisible(dev.off())


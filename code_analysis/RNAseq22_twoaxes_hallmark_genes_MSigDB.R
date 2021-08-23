
# Explore gene sets first 

msigdbr_species()

msigdbr_collections() %>% as.data.frame()

h_gene_sets = msigdbr(species = "Mus musculus", category = "H")

head(h_gene_sets)

h_gene_sets$gs_name %>% unique()


h_gene_sets %>% 
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") # 200 genes 
# but remember, these are genes associated with LPS challenge or similiar stimulation. Does it reflect 'chronic' inflammation level? 

h_gene_sets %>% 
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") # 199 genes 

h_gene_sets %>% 
  filter(grepl("Cyp",gene_symbol)) %>% 
  .$gs_name %>% 
  unique()

# Immune hallmark gene (C7)
# how do I know C7 is immune hallmark gene set? 
# Go to; http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C7
C7_immune_gene_sets = msigdbr(species = "Mus musculus", category = "C7")
head(C7_immune_gene_sets)

C7_immune_gene_sets$gs_name %>% unique() # realizing there are tooo many subcat


#  of most interest =============================================================


my_tissue = "Liver"
my_tissue = "Spleen"
if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}_second.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}_second.RDS"))}


hallmark = "HALLMARK_INFLAMMATORY_RESPONSE"
my_hallmark = "Inflammatory Response" # for title 

hallmark = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
my_hallmark = "Oxidative phosphorylation"

h_gene_sets %>% 
  filter(gs_name == hallmark) %>% 
  .$gene_symbol -> selected_genes


limma_list$status %>% 
  filter(symbol %in% selected_genes) %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_status = logFC,
         pval_status = P.Value) -> status_df

limma_list$cort %>% 
  filter(symbol %in% selected_genes) %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_cort = logFC,
         pval_cort = P.Value) -> cort_df

df <- full_join(status_df, cort_df) %>% 
  mutate(Sig = ifelse(pval_status >= 0.05 & pval_cort >= 0.05, "None", 
                      ifelse(pval_status < 0.05 & pval_cort < 0.05,"Both",
                             ifelse(pval_status < 0.05,"Status-specific","CORT-specific")))) %>% 
  mutate(Sig = factor(Sig, levels =c("Both","Status-specific","CORT-specific","None")))

max(abs(c(df$logFC_status, df$logFC_cort)))*1.05 -> lim

png(filename = glue("results_figures/twoaxis_{hallmark}_{my_tissue}second.png"),
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
  geom_text_repel(data = df %>% 
               filter(Sig == "Both"), 
             aes(logFC_status, logFC_cort,label = symbol),
            color ="black",
             size = 2.5)+
  
    labs(x = "Status effect (Sub <-> Dom)",
       y = "CORT effect (Low <-> High)",
       color = "",
       title = glue("{my_tissue}: Hallmark {my_hallmark}"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "top")+
  scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))

invisible(dev.off())


#  Now automate it =============================================================


my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}


h_gene_sets$gs_name %>% unique() -> hm_names
lim = 1.5

for(i in 1:length(hm_names)){
    my_hallmark <- gsub("HALLMARK_","",hm_names[i])
    hallmark <- hm_names[i]
  h_gene_sets %>% 
    filter(gs_name == hallmark) %>% 
    .$gene_symbol -> selected_genes
  
  limma_list$status %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    rename(logFC_status = logFC,
           pval_status = P.Value) -> status_df
  
  limma_list$cort %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    rename(logFC_cort = logFC,
           pval_cort = P.Value) -> cort_df
  
  df <- full_join(status_df, cort_df) %>% 
    mutate(Sig = ifelse(pval_status >= 0.05 & pval_cort >= 0.05, "None", 
                        ifelse(pval_status < 0.05 & pval_cort < 0.05,"Both",
                               ifelse(pval_status < 0.05,"Status-specific","CORT-specific")))) %>% 
    mutate(Sig = factor(Sig, levels =c("Both","Status-specific","CORT-specific","None")))
  
  
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
         title = glue("{my_tissue}: {my_hallmark}"))+
    theme_bw(base_size = 7)+
    theme(legend.position = "top")+
    scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
    scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0))) -> temp
  
  png(filename = glue("results_figures/hallmark_geneset_plots/twoaxis_{hallmark}_{my_tissue}.png"),
      width = 7, height = 7.6, units = "cm", res = 600)
  print(temp)
  dev.off()
  
}


# other pathways ==============================================================================
msigdbr_collections() %>% as.data.frame()

# -----------------------------------------------------
msigdbr(species = "Mus musculus", category = "C7") %>% 
  filter(grepl("CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE",gs_name,ignore.case = T)) -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory ="GO:BP") %>% 
  filter(grepl("oxidative",gs_name,ignore.case = T)) -> my_gs_sets


msigdbr(species = "Mus musculus", subcategory ="GO:BP") %>% 
  filter(grepl("fatty_acid_biosyn",gs_name,ignore.case = T)) -> my_gs_sets
# -------------------------------------------------------------

# immu
# inflamm
my_gs_sets %>% 
  .$gs_name %>% unique() -> my_gs_names
my_gs_names
my_gs_names <- my_gs_names[c(1)]

my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}



for(i in 1:length(my_gs_names)){
  my_hallmark <- gsub("GO_","",my_gs_names[i])
  hallmark <- my_gs_names[i]
  my_gs_sets %>% 
    filter(gs_name == hallmark) %>% 
    .$gene_symbol -> selected_genes
  
  limma_list$status %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    rename(logFC_status = logFC,
           pval_status = P.Value) -> status_df
  
  limma_list$cort %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    rename(logFC_cort = logFC,
           pval_cort = P.Value) -> cort_df
  
  df <- full_join(status_df, cort_df) %>% 
    mutate(Sig = ifelse(pval_status >= 0.05 & pval_cort >= 0.05, "None", 
                        ifelse(pval_status < 0.05 & pval_cort < 0.05,"Both",
                               ifelse(pval_status < 0.05,"Status-specific","CORT-specific")))) %>% 
    mutate(Sig = factor(Sig, levels =c("Both","Status-specific","CORT-specific","None")))
  
  max(abs(c(df$logFC_status, df$logFC_cort)))*1.05 -> lim
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
    geom_text_repel(data = df %>% 
                      filter(Sig == "Status-specific"), 
                    aes(logFC_status, logFC_cort,label = symbol),
                    color ="black",
                    size = 2.5)+
    labs(x = "Status effect (Sub <-> Dom)",
         y = "CORT effect (Low <-> High)",
         color = "",
         title = glue("{my_tissue}: {my_hallmark}"))+
    theme_bw(base_size = 7)+
    theme(legend.position = "top")+
    scale_x_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0)))+
    scale_y_continuous(limits = c(-lim,lim), expand=expansion(mult=c(0,0.0))) -> temp
  
  png(filename = glue("results_figures/GO_geneset_plots/twoaxis_{hallmark}_{my_tissue}.png"),
      width = 7, height = 7.6, units = "cm", res = 600)
  print(temp)
  dev.off()
  
}

my_tissue = "Liver"
my_tissue = "Spleen"
if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}

my_logFC_threshold = 0.15


genecat = "Hallmark"
hallmark_sets <- readRDS("results_RNAseqRDS/hallmark_sets.RDS")
hallmark_sets %>% 
  filter(grepl("OXIDATIVE",gs_name,ignore.case = T)) -> my_gs_sets

my_gs_sets %>% 
  .$gs_name %>% unique() -> my_gs_names

my_gs_names
my_titles <- c("Oxidative phosphorylation")

i = 1

my_hallmark <- my_titles[i]
hallmark <- my_gs_names[i]
my_gs_sets %>% 
  filter(gs_name == hallmark) %>% 
  .$gene_symbol -> selected_genes



msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("DNA_repair",gs_name,ignore.case = T)) -> my_gs_sets
my_gs_sets %>% 
  .$gs_name %>% unique() -> my_gs_names
my_gs_names

my_gs_names <- my_gs_names[c(3)]
my_titles <- c("DNA repair")

my_hallmark <- my_titles[i]
hallmark <- my_gs_names[i]
my_gs_sets %>% 
  filter(gs_name == hallmark) %>% 
  .$gene_symbol -> selected_genes
# my_gs_sets %>% 
#   filter(gs_name == hallmark) %>% 
#   .$gene_symbol -> selected_genes_repair





limma_list$cort %>% # CORT !!! 
  filter(symbol %in% selected_genes) %>% 
  select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC > my_logFC_threshold,"Positive association","Negative association"))) %>%
  mutate(Sig = factor(Sig, levels = c("Positive association","Negative association", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df


keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'grey'))

table(keyvals) %>% 
  as.data.frame() %>% 
  mutate(prop = glue("{round(Freq/sum(Freq),3)*100}")) %>% 
  mutate(text = glue("{Freq} genes ({prop}%)")) -> my_texts

if(is.na(my_texts[3,4])){my_texts[3,4] <- "0 genes"}


keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'


df %>% 
  filter(P.Value <0.05) %>%
  filter(abs(logFC) > my_logFC_threshold) %>%
  arrange(logFC) %>%
  mutate(effects = abs(sin(logFC)*-log10(P.Value))) %>%
  top_frac(.,0.3) %>%
  # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
  .$symbol -> for_label


EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue}: {my_hallmark}"),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                # drawConnectors = TRUE, # just for a few plots 
                # widthConnectors = 0.5, # just for a few plots 
                vline = c(-0.15, 0.15),
                vlineCol = c('grey90'),
                vlineType = c( 'dashed'),
                vlineWidth = c(0.3),
                hline = c(0.05),
                hlineCol = c('grey90'),
                hlineType = c( 'dashed'),
                hlineWidth = c(0.3),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 2.1)+
  annotate("text", x =0.8, y = 4.4,color = "purple4" , size = 2.2,
           label = glue("Positively associated with CORT \n{my_texts[3,4]}"))+
  annotate("text", x = -0.8, y = 4.4,color = "orange" , size = 2.2,
           label = glue("Negatively associated with CORT \n{my_texts[2,4]}"))+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
  theme_bw(base_size = 7)+
  labs(color = "",
       caption = paste0('total = ', nrow(df), ' genes'),
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.title = element_text(size = 11),
        plot.subtitle = element_blank()) -> temp



png(filename = glue("manuscript/figures/fig3/EnVol/by_CORT_{my_tissue}_{my_hallmark}.png"),
    width = 8, height = 7.3, units = "cm", res = 600)

print(temp)
invisible(dev.off())



df %>% filter(Sig == "Positive association")


# IDEA = ======================================================================================
# I am borrowing moduleEigengenes(datExpr, moduleColors) function to calculate Eigengenes of selected gene sets 

my_networkType = "signed hybrid"
my_tissue = "Liver"
my_power= 5
my_nrow = 2 # for liver 
my_height = 5 # for liver 



my_networkType = "signed hybrid"
my_tissue = "Spleen"
my_power= 14
my_nrow = 1 # for spleen
my_height = 2.5 # for spleen 


datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
dim(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))


rownames(datExpr)
colnames(alldata)
datTraits <- rnaseq_sampleid %>% 
  filter(sampleID %in% rownames(datExpr)) %>% 
  left_join(alldata) %>% 
  mutate(status = ifelse(domgroup == "Alpha", 1,0)) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", 1,0)) %>% 
  column_to_rownames("sampleID") %>% 
  dplyr::select(
    GD01_BW, GD14_bw, spleen_weight_mg, AGD, pair_status,
    status,ds, glicko_rank, win_prop, loss_prop, ind_cert,
    preem_flee_perc, flee_perc, sub_pos_perc, freeze_perc, 
    despotism, cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
    `B cells_GD14`, `CD4_CD8_ratio_GD14`, `Cytotoxic T_GD14`, `DCs_GD14`,
    `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
    `Neutro_GD14`, `NK cells_GD14`, `NLR_GD14`, `T cells_GD14`
  )


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# now my code 
grcm38 %>% 
  filter(symbol %in% selected_genes) %>% 
  .$ensgene -> selected_ensgene

colnames(datExpr) %in% selected_ensgene -> temp_vec
for(i in 1:length(temp_vec)){
  if(temp_vec[i] == T){temp_vec[i] <- "Oxi"} else{temp_vec[i] <-"Out"}  
}
temp_vec
grcm38 %>% 
  filter(symbol %in% selected_genes_repair) %>% 
  .$ensgene -> selected_ensgene_repair

colnames(datExpr) %in% selected_ensgene_repair -> temp_vec_repair
for(i in 1:length(temp_vec_repair)){
  if(temp_vec[i] == "Oxi"){temp_vec[i] <- "Oxi"} else{
    if(temp_vec_repair[i] == T){temp_vec[i] <- "DNArep"} else{temp_vec[i] <-"Out"}    
  }
  
}



temp_vec

In_and_Out = moduleEigengenes(datExpr,temp_vec)$eigengenes


cbind(MEs,In_and_Out) %>%
  rownames_to_column("sampleID") %>%
  left_join(datTraits %>%
              rownames_to_column("sampleID")) %>% 
  mutate(status = ifelse(status ==1, "Alpha", "Subordinate"))-> ME_df

ME_df %>% 
  ggplot(aes(`T cells_GD14`, MEIn, color = status))+
  geom_point()+
  geom_smooth(method = "lm")

ME_df %>% 
  ggplot(aes(MEbrown, MEIn,color = status, fill = status))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm")

ME_df %>% 
  ggplot(aes(cort_post, MEturquoise,color = status, fill = status))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm")

ME_df %>% 
  ggplot(aes(MEbrown, MEgreen,color = status, fill = status))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm")


ME_df %>% 
  left_join(rnaseq_sampleid) %>% 
  ggplot(aes(cort_post, MEblue,color = status))+
  geom_point()+
  geom_text(aes(label = subjectID))+
  geom_smooth(method = "lm")
rnaseq_sampleid

summary(lm(MEgreen~MEIn, data = ME_df))



ME_df %>% 
  ggplot(aes(cort_post, `B cells_GD14`))+
  geom_point()+
  geom_smooth(method = "lm")




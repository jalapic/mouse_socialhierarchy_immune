# install.packages("msigdbr")

# Explore gene sets first 
library(msigdbr)

msigdbr_collections() %>% as.data.frame()


# Immune hallmark gene (C7)
# how do I know C7 is immune hallmark gene set? 
# Go to; http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C7

my_logFC_threshold = 0.15
# other pathways ==============================================================================
msigdbr_collections() %>% as.data.frame()

# msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
# saveRDS(reactome_sets,"results_RNAseqRDS/reactome_sets.RDS")
genecat = "Reactome"
reactome_sets <- readRDS("results_RNAseqRDS/reactome_sets.RDS")
reactome_sets %>% 
  filter(grepl("immune_system",gs_name,ignore.case = T)) -> my_gs_sets
reactome_sets %>% 
  filter(grepl("oxidative",gs_name,ignore.case = T)) -> my_gs_sets
reactome_sets %>% 
  filter(grepl("metabolism",gs_name,ignore.case = T)) -> my_gs_sets

genecat = "Hallmark"

# msigdbr(species = "Mus musculus", category ="H") -> hallmark_sets
# saveRDS(hallmark_sets,"results_RNAseqRDS/hallmark_sets.RDS")
hallmark_sets <- readRDS("results_RNAseqRDS/hallmark_sets.RDS")
hallmark_sets %>% 
  filter(grepl("oxidative",gs_name,ignore.case = T)) -> my_gs_sets
hallmark_sets %>% 
  filter(grepl("inflam",gs_name,ignore.case = T)) -> my_gs_sets

genecat = "GO"
msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("fatty_acid",gs_name,ignore.case = T)) -> my_gs_sets


my_gs_sets %>% 
  .$gs_name %>% unique() -> my_gs_names
my_gs_names
my_gs_names <- my_gs_names[c(2,3,6,7)]
my_gs_names

my_tissue = "Liver"
my_tissue = "Spleen"

limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
limma_vdl <- readRDS(glue("results_RNAseqRDS/limma_vdl_{my_tissue}.RDS"))


for(i in 1:length(my_gs_names)){

my_hallmark <- gsub(glue("{genecat}_"),"",my_gs_names[i],ignore.case = T)
hallmark <- my_gs_names[i]
my_gs_sets %>% 
  filter(gs_name == hallmark) %>% 
  .$gene_symbol -> selected_genes

limma_list$status %>% 
  filter(symbol %in% selected_genes) %>% 
  select(symbol, logFC, P.Value) %>% 
  mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                      ifelse(logFC>my_logFC_threshold,"Alpha genes","Subordinate genes"))) %>%
  mutate(Sig = factor(Sig, levels = c("Alpha genes","Subordinate genes", "N.S."))) %>% 
  mutate(P.Value = ifelse(P.Value == 0, 1/10000,P.Value)) %>% 
  unique() -> df


keyvals <- ifelse(
  df$logFC < -my_logFC_threshold & df$P.Value<0.05, 'orange',
  ifelse(df$logFC > my_logFC_threshold & df$P.Value<0.05, 'purple4',
         'black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'
table(keyvals)

png(filename = glue("results_figures/geneset_volcano/EnVol_{my_tissue}_{my_hallmark}.png"),
    width = 12, height = 11, units = "cm", res = 600)
df %>% 
  filter(P.Value <0.05 & abs(logFC) >0.25) %>% 
  .$symbol -> for_label
EnhancedVolcano(df,
                lab = df$symbol,
                selectLab= for_label,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue} DEG by Status: {genecat}-{my_hallmark}"),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "black",
                labSize = 2)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.5),breaks = c(0,1,2,3,4))+
  theme_bw(base_size = 7)+
  labs(color = "",
       y = bquote(~-Log[10]~italic(eFDR)))+
  theme(legend.position = "none",
        plot.subtitle = element_blank()) -> temp
print(temp)
invisible(dev.off())
}

# for manuscript  ====================================================================



# 1 . GO term  - immune 
genecat = "GO"
msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("immune",gs_name,ignore.case = T)) -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("monocarboxyl",gs_name,ignore.case = T)) %>% 
  filter(grepl("process",gs_name,ignore.case = T))   -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("ic_process",gs_name,ignore.case = T)) %>%
  filter(grepl("fatty_acid",gs_name,ignore.case = T))   -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("stress",gs_name,ignore.case = T)) %>% 
  filter(grepl("oxidative",gs_name,ignore.case = T))   -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  # filter(grepl("stress",gs_name,ignore.case = T)) %>% 
  filter(grepl("oxidoreductase",gs_name,ignore.case = T))   -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("DNA_repair",gs_name,ignore.case = T)) -> my_gs_sets

msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("inflammat",gs_name,ignore.case = T)) -> my_gs_sets

genecat = "GO"
msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("response_to_oxidative_stress",gs_name,ignore.case = T)) -> my_gs_sets


genecat = "GO"
msigdbr(species = "Mus musculus", subcategory = "GO:BP") %>% 
  filter(grepl("mitotic_cell_cycle",gs_name,ignore.case = T)) -> my_gs_sets







my_gs_sets %>% 
  .$gs_name %>% unique() -> my_gs_names
my_gs_names








# my_gs_names <- my_gs_names[c(3,12,21,78,79,81)]
my_gs_names <- my_gs_names[c(3,22)]
my_titles <- c("Innate immune response","Adaptive immune response")




my_titles <- c("Monocarboxylic acid biosythesis","Monocarboxylic acid catabolism","Monocarboxylic acid metabolism")

my_gs_names <- my_gs_names[c(1,2,18)]
my_titles <- c("Fatty acid biosynthesis","Fatty acid catabolism", "Reg. of fatty acid biosynthesis")

my_gs_names <- my_gs_names[c(18)]
my_titles <- c("Response to oxidative stress")

my_gs_names <- my_gs_names[c(3)]
my_titles <- c("Oxidoreductase activity")

# my_gs_names <- my_gs_names[c(3,7,8,9)]
# my_titles <- c("DNA repair","Negative reg. of DNA repair", 
#                "Positive reg. of DNA repair","Reg of DNA repair")
my_gs_names <- my_gs_names[c(3)]
my_titles <- c("DNA repair")



my_gs_names <- my_gs_names[c(8,10)]
my_titles <- c("Reg. of response to oxidative stress", 
               "Response to oxidative stress")


my_gs_names <- my_gs_names[c(1,8)]
my_titles <- c("Mitotic cell cycle", 
               "Reg. of mitotic cell cycle")


# 2. Hallmark - inflammatory and oxidative stress
genecat = "Hallmark"
# msigdbr(species = "Mus musculus", category ="H") -> hallmark_sets
# saveRDS(hallmark_sets,"results_RNAseqRDS/hallmark_sets.RDS")
hallmark_sets <- readRDS("results_RNAseqRDS/hallmark_sets.RDS")

hallmark_sets %>%
  filter(grepl("oxidative",gs_name,ignore.case = T)|grepl("glycolysis",gs_name,ignore.case = T)) -> my_gs_sets
# hallmark_sets %>%
#   filter(grepl("tgf",gs_name,ignore.case = T)) -> my_gs_sets


my_gs_sets %>%
  .$gs_name %>% unique() -> my_gs_names
my_gs_names


my_titles <- c("Glycolysis","Oxidative phosphorylation")

# my_titles <- c("Tgf-beta signaling")


# 3. Reactome
genecat = "Reactome"
reactome_sets <- readRDS("results_RNAseqRDS/reactome_sets.RDS")
reactome_sets %>%
  filter(grepl("immune_system",gs_name,ignore.case = T)) -> my_gs_sets

my_gs_sets %>%
  .$gs_name %>% unique() -> my_gs_names
my_gs_names

my_gs_names <- my_gs_names[c(1,2,4)]
my_titles <- c("Adaptive immune system", "Cytokine signaling","Innate immune system")


# # 4. imprinted gene (also do for brain later)
# msigdbr(species = "Mus musculus", category ="C2") %>%
#   filter(grepl("bri",gs_name,ignore.case = T)) %>%
#   filter(grepl("imprinted",gs_name,ignore.case = T))-> my_gs_sets
# 
# my_gs_sets %>%
#   .$gs_name %>% unique() -> my_gs_names
# 
# 
# my_gs_sets$gs_name %>% unique()
# 
# 
# my_titles <- "Imprinted genes"

# Plot it ====================================================================
my_logFC_threshold = 0.15

my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}



for(i in 1:length(my_gs_names)){
  
  my_hallmark <- my_titles[i]
  hallmark <- my_gs_names[i]
  my_gs_sets %>% 
    filter(gs_name == hallmark) %>% 
    .$gene_symbol -> selected_genes
  
  limma_list$status %>% 
    filter(symbol %in% selected_genes) %>% 
    select(symbol, logFC, P.Value) %>% 
    mutate(Sig = ifelse(P.Value >= 0.05, "N.S.",
                        ifelse(logFC>my_logFC_threshold,"Alpha genes","Subordinate genes"))) %>%
    mutate(Sig = factor(Sig, levels = c("Alpha genes","Subordinate genes", "N.S."))) %>% 
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

  if(is.na(my_texts[2,4])){my_texts[2,4] <- "0 genes"}
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
    top_n(40) %>% 
    # top_frac(.,1) %>%
    # filter(P.Value <0.01 | abs(logFC) >0.5) %>%
    .$symbol -> for_label
  
  png(filename = glue("manuscript/figures/fig3/EnVol/{my_tissue}_{my_hallmark}_biggerfont.png"),
      width = 8, height = 7.3, units = "cm", res = 600)
 
  
  
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
                  # widthConnectors = 0.1, # just for a few plots
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
                  pointSize = 2.6,
                  labCol = "black",
                  labSize = 2.4)+
    annotate("text", x = 1, y = 4.4,color = "purple4" , size = 3,
             label = glue("Upregulated in alpha \n{my_texts[3,4]}"))+
    annotate("text", x = -0.75, y = 4.4,color = "orange" , size = 3,
             label = glue("Upregulated in subordinate \n{my_texts[2,4]}"))+
    scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
    scale_y_continuous(limits = c(-0.1,4.8),breaks = c(0,1,2,3,4),expand=expansion(mult=c(0.0005,0.0)))+
    theme_bw(base_size = 7)+
    labs(color = "",
         caption = paste0('total = ', nrow(df), ' genes'),
         y = bquote(~-Log[10]~italic(eFDR)))+
    theme(legend.position = "none",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 9),
          plot.caption = element_text(size = 9),
          plot.subtitle = element_blank()) -> temp
  print(temp)
  invisible(dev.off())
}


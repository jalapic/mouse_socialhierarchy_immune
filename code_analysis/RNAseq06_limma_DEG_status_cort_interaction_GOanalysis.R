my_tissue = "Spleen"
my_tissue = "Liver"
my_ont = "BP"

my_showCategory = 100
my_logFC_threshold = 0.15
limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS")) %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

#colnames(limma_list$status)
limma_list %>% 
  map(~summarise(.,Up = sum(logFC>0),
                 Down = sum(logFC<0))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_status_DEG <- limma_list$status 
limma_cort_DEG <- limma_list$cort 
# limma_interaction_DEG <- limma_list$interaction 

# limma_list %>% map(~hist(.$logFC))

gettop10GO <- function(limma_df,my_showCategory = 10){
  go_df_up <- limma_df %>% 
    filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  
  go_df_down <- limma_df %>% 
    filter(logFC<0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  ggo_up <- enrichGO(gene = go_df_up$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.2,
                     qvalueCutoff  = 0.20)
  
  ggo_down <- enrichGO(gene = go_df_down$entrez %>% unique() ,
                       OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = my_ont,
                       readable = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.2, # change to 0.2 for Liver-status only
                       qvalueCutoff  = 0.20) 
  
  fortify(
    ggo_up,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Up") -> temp1
  
  
  fortify(
    ggo_down,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Down") -> temp2
  return(rbind(temp1,temp2))
  
}


gettop10GO(limma_status_DEG, my_showCategory) %>% 
  mutate(comparison = "Social status") -> top10go1

gettop10GO(limma_cort_DEG, my_showCategory ) %>% 
  mutate(comparison = "Corticosterone level (GD14)") -> top10go2

# gettop10GO_interaction(limma_interaction_DEG, my_showCategory ) %>% 
#   mutate(comparison = "Interaction (CORT * Status)") -> top10go3

rbind(top10go1,top10go2) %>% 
  mutate(tissue = my_tissue) -> top10_GOterms

write.csv(top10_GOterms,glue("results_tables/top{my_showCategory}_GOterms_abovelogFC{my_logFC_threshold}_{my_tissue}.csv"), row.names = F)









View(top10_GOterms %>% select(Description, direction, comparison))



saveRDS(top10_GOterms,glue("results_RNAseqRDS/top{my_showCategory}_GOterms_abovelogFC{my_logFC_threshold}_{my_tissue}_newnew.RDS"))


# view

top10_GOterms$Description[1:5]
top10_GOterms$Description[6:10]
top10_GOterms$Description[11:15]
top10_GOterms$Description[16:20]
top10_GOterms$Description[21:25]


limma_status_DEG %>% 
  arrange(desc(logFC)) %>% 
  head(5) %>% .$symbol

limma_status_DEG %>% 
  arrange(logFC) %>% 
  head(20) %>% .$symbol

limma_cort_DEG%>% 
  arrange(desc(logFC)) %>% 
  head(20) %>% .$symbol

limma_cort_DEG %>% 
  arrange(logFC) %>% 
  head(10) %>% .$symbol

limma_interaction_DEG %>% 
  arrange(desc(abs(logFC))) %>% 
  head(6) %>% .$symbol




gettop10GO_interaction <- function(limma_df,my_showCategory = 10){
  go_df <- limma_df %>% 
    # filter(logFC>0) %>% 
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez))
  
  ggo <- enrichGO(gene = go_df$entrez %>% unique(),
                  OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = my_ont,
                  readable = T,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.20)
  
  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(direction = "Both") -> temp1
  
  return(temp1)
  
}
# Plot them =============================================================================
png(glue("results_figures/enrichGO_{my_tissue}_limma_upregulated{my_compare2}_{my_ont}.png"),
    width = 9, height = 4, units = 'in', res = 300)

clusterProfiler::dotplot(ggo_up, orderBy = "Count")+
  ggtitle(glue("{my_tissue}: Upregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()


dev.off()

png(glue("results_figures/enrichGO_{my_tissue}_limma_downregulated{my_compare2}_{my_ont}.png"), 
    width = 9, height = 4, units = 'in', res = 300)

clusterProfiler::dotplot(ggo_down, x = 'GeneRatio')+
  ggtitle(glue("{my_tissue}: Downregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()

dev.off()

# Extract info to make a summary table ===============================================
length(ggo_up@gene) # upregulated genes 
length(ggo_down@gene) # downregulated genes 



ggo_up %>%
  arrange(desc(Count),p.adjust) %>%
  mutate(p.adjust = round(p.adjust,3)) %>%
  mutate(qvalue = round(qvalue,3)) %>%
  as.data.frame() %>% 
  select(Description,Count,p.adjust,geneID) %>% 
  mutate(direction = "Up") -> temp1

ggo_down %>%
  arrange(desc(Count),p.adjust) %>%
  mutate(p.adjust = round(p.adjust,3)) %>%
  mutate(qvalue = round(qvalue,3)) %>% 
  as.data.frame() %>% 
  select(Description,Count,p.adjust,geneID) %>% 
  mutate(direction = "Down")-> temp2

rbind(temp1,temp2) %>% 
  mutate(tissue = my_tissue) %>% 
  mutate(comparison = my_compare2) -> temp
summary_GOterms %>% rbind(temp) -> summary_GOterms


tail(summary_GOterms)

saveRDS(summary_GOterms, glue("results_RNAseqRDS/summary_GOterms_{my_tissue}.RDS"))

# Interaction - rather all together? 

summary_GOterms %>% 
  filter(direction == "Down") %>% 
  filter(comparison == "Interaction")


fortify(
  ggo_up,
  showCategory = 10,
  by = "Count",
  split = NULL,
  includeAll = TRUE
) %>% 
  arrange(desc(GeneRatio)) %>% 
  .$Description

# KEGG gene set enrichment analysis ================================================

search_kegg_organism('Mus', by='scientific_name')
search_kegg_organism('mmu', by='kegg_code')



kk_up <- enrichKEGG(gene    = go_df_up$entrez,
                      organism     = 'mmu',
                      pvalueCutoff = 0.05)
head(kk_up)

clusterProfiler::dotplot(kk_up, orderBy = "Count")+
  ggtitle(glue("KEGG - {my_tissue}: Upregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()


kk_down <- enrichKEGG(gene    = go_df_down$entrez,
                      organism     = 'mmu',
                      pvalueCutoff = 0.05)
head(kk_down)

clusterProfiler::dotplot(kk_down, x = 'GeneRatio')+
  ggtitle(glue("{my_tissue}: Downregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()


d <- go_df_up %>% as.data.frame()
geneList <- d[,2]
names(geneList) <- d[,1]
geneList <- sort(geneList, decreasing = TRUE)

kk2_up <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               nPerm        = 1000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid",
               verbose = F)

png(glue("results_figures/gseKEGG_{my_tissue}_limma_upregulated{my_compare2}.png"), 
    width = 9, height = 4, units = 'in', res = 300)

clusterProfiler::dotplot(kk2_up, orderBy = "GeneRatio")+
  ggtitle(glue("KEGG - {my_tissue}: Upregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()

dev.off()

d <- go_df_down %>% as.data.frame()
geneList <- d[,2]
names(geneList) <- d[,1]
geneList <- sort(geneList, decreasing = TRUE)

kk2_down <- gseKEGG(geneList     = geneList,
                  organism     = 'mmu',
                  nPerm        = 1000,
                  minGSSize    = 3,
                  maxGSSize    = 800,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "none",
                  keyType       = "ncbi-geneid",
                  verbose = F)

png(glue("results_figures/gseKEGG_{my_tissue}_limma_downregulated{my_compare2}.png"), 
    width = 9, height = 4, units = 'in', res = 300)
clusterProfiler::dotplot(kk2_down, orderBy = "GeneRatio")+
  ggtitle(glue("KEGG - {my_tissue}: Downregulated (~{my_compare})"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()
dev.off()


head(kk2_down,10)

kk2_down %>% as.data.frame() %>%
  filter(Description == "Th1 and Th2 cell differentiation") %>%
  .$core_enrichment %>%
  str_split(.,"/") %>%
  as.data.frame(col.names = "entrez") %>%
  mutate(entrez = as.integer(entrez)) %>%
  left_join(grcm38 %>% select(entrez, symbol, description))



# For genewalk ========================================================

my_tissue = "Spleen"
my_tissue = "Liver"

my_logFC_threshold = 0.0

limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS")) %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) %>% 
  map(~select(.,entrez,logFC)) 

limma_list$status %>% 
  filter(logFC>0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_status_UP.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)


limma_list$status %>% 
  filter(logFC<0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_status_DOWN.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)


limma_list$cort %>% 
  filter(logFC>0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_cort_UP.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)


limma_list$cort %>% 
  filter(logFC<0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_cort_DOWN.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)



limma_list$interaction %>% 
  filter(logFC>0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_interaction_UP.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)


limma_list$interaction %>% 
  filter(logFC<0) %>% 
  select(entrez) %>% 
  write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_interaction_DOWN.csv"),
              row.names = F, quote = F, sep = "\t", col.names = F)

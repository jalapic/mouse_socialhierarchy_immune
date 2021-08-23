# run step 2 first! 


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

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

modNames = substring(names(MEs), 3)

# module = "pink"
# module = "blue"
# module = "cyan"
# module = "grey60"
#
# module_list = c("pink","cyan", "blue", "grey60") # For liver

# module = "green"
# module = "brown"
#
# module_list = c("green", "brown") # For spleen ~ with cort GD14


## ==============================================
gettop10GO_WGCNA <- function(module,my_showCategory = 10){
  
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    select(entrez) -> go_df_wgcna
  
  
  ggo <- enrichGO(gene = go_df_wgcna$entrez %>% unique(),
                     OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = my_ont,
                     readable = T,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.50)
  

  fortify(
    ggo,
    showCategory = my_showCategory,
    by = "Count",
    split = NULL,
    includeAll = TRUE
  ) %>% 
    arrange(desc(GeneRatio)) %>% 
    mutate(module = module) -> temp1

  return(rbind(temp1))
  
}

my_ont = "BP"
my_showCategory = 100

# gettop10GO_WGCNA("pink")
# gettop10GO_WGCNA("cyan")


WGCNA_GOs <- vector('list', length(allcolors))

moduleColors %>% unique() -> allcolors

for(i in 1:length(allcolors)){
  gettop10GO_WGCNA(allcolors[i],my_showCategory) -> WGCNA_GOs[[i]]
}

WGCNA_GOs %>% 
  rbindlist() -> wgcna_all_gos

write.csv(wgcna_all_gos, 
          glue("results_tables/wgcna_all_gos_{my_tissue}_{my_showCategory}.csv"),
          row.names = F)





saveRDS(wgcna_all_gos,glue("results_RNAseqRDS/wgcna_all_gos_{my_tissue}.RDS"))


# my_tissue = "Liver"
# my_tissue = "Spleen"
# wgcna_all_gos <- readRDS(glue("results_RNAseqRDS/wgcna_all_gos_{my_tissue}.RDS"))
# View(wgcna_all_gos %>% select(module, Description))







# GO over-representation analysis ================================================

ggo_wgcna <- enrichGO(gene = go_df_wgcna$entrez ,
                   OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = my_ont,
                   readable = T,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.10)



png(glue("results_figures/enrichGO_{my_tissue}_WGCNA_{module}_{my_ont}.png"),
    width = 9, height = 4, units = 'in', res = 300)

clusterProfiler::dotplot(ggo_wgcna, orderBy = "Count")+
  ggtitle(glue("{my_tissue}: {module} Module"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()


dev.off()

# KEGG over-representation enrichment analysis ================================================

kk_wgcna <- enrichKEGG(gene    = go_df_wgcna$entrez,
                    organism     = 'mmu',
                    pvalueCutoff = 0.05)
head(kk_wgcna)

png(glue("results_figures/gseKEGG_{my_tissue}_WGCNA_{module}.png"), 
    width = 9, height = 4, units = 'in', res = 300)
clusterProfiler::dotplot(kk_wgcna)+
  ggtitle(glue("KEGG - {my_tissue}: Module {module}"))+
  theme(plot.title = element_text(size = 12))+
  scale_color_viridis()

dev.off()

# entrez list for genewalk =========================================


moduleColors %>% unique() -> allcolors
for(i in 1:length(allcolors)){
  allcolors[i] -> module
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    filter(ensgene %in% module_gene) %>% 
    filter(!is.na(entrez)) %>% 
    select(entrez) %>% 
    write.table(glue("results_GeneWalk/forGeneWalk_{my_tissue}_WGCNA_{module}.csv"), sep = "\t",
                quote = F,
                col.names = F,
                row.names = F)
}


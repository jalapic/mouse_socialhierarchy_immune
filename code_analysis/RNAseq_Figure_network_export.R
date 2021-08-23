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

datExpr <- readRDS(glue("C:/Users/PSYC-wl8856/Desktop/Github/immune/results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
net <- readRDS(glue("C:/Users/PSYC-wl8856/Desktop/Github/immune/results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))

# select module
my_module = "pink"
my_module = "cyan"
my_module = "blue"

my_module = "brown"
my_module = "green"

# ====================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleColors %>% table()
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
modNames = substring(names(MEs), 3)

modules = my_module

# Select module probes
probes = colnames(datExpr)

inModule = is.finite(match(moduleColors, modules));

modProbes = probes[inModule];

modGenes = grcm38$symbol[match(modProbes, grcm38$ensgene)];
length(modGenes)

cbind(MODULE = my_module, GENE = modGenes) %>% 
  as.data.frame() %>% 
  filter(!is.na(GENE)) %>% 
  mutate(OVERLAP = MODULE,
         NODE = GENE) -> genesets

selected <- (moduleColors == my_module)
TOM = TOMsimilarityFromExpr(datExpr[,selected],
                            power = my_power,
                            networkType = "signed hybrid")

dimnames(TOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])
cyt$edgeData %>% 
  mutate(HEAD = fromNode,
         TAIL = toNode,
         WEIGHT = weight) %>% 
  select( HEAD,TAIL, WEIGHT) %>%
  filter(!is.na(TAIL)) %>% 
  filter(!is.na(HEAD)) -> ass

write.table(ass,glue("results_networks/WGCNA{my_module}_{my_tissue}.txt"),
            row.names = F, quote = F, sep = "\t")
write.table(modGenes,glue("results_networks/genes_{my_module}_{my_tissue}.txt"),
            row.names = F, quote = F, sep = "\t", col.names = F)

# WGCNA - intramodular connectivity ==========================

ADJ=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ, moduleColors) 

Alldegrees1 %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) -> kIN_df

xx <- cbind(ensgene = colnames(datExpr), module = moduleColors) %>% 
  as.tibble()


kIN_df %>% 
  left_join(xx) -> kIN_df

saveRDS(kIN_df,glue("results_RNAseqRDS/kIN_dataframe_{my_tissue}.RDS"))



# =========================================================
if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
limma_result_of_interest <- limma_list$status} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
limma_result_of_interest <- limma_list_spleen$cort }


limma_result_of_interest %>% 
  select(symbol,logFC) %>% 
  write.table(.,glue("results_networks/for_ppi_liver_status.txt"),
              row.names = F, quote = F, sep = "\t", col.names = F)
  

limma_list_spleen$cort %>% 
  select(symbol,logFC) %>% 
  write.table(.,glue("results_networks/for_ppi_spleen_cort.txt"),
              row.names = F, quote = F, sep = "\t", col.names = F)

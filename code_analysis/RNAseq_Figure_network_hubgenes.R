
overlap_list <- list()
both_list <- list()
bothbutoverlap_list <- list()




my_tissue = "Liver"
my_module = "pink"
cutoff = round(176*0.15)

my_tissue = "Liver"
my_module = "cyan"
cutoff = round(71*0.15)
# name_any_gene("Acly")
# name_any_gene("Pemt")
# name_any_gene("Srebf1")
# name_any_gene("Prdx4") # interesting!! https://www.mdpi.com/1422-0067/19/9/2509

my_tissue = "Liver"
my_module = "blue"
cutoff = round(730*0.05)





my_tissue = "Spleen"
my_module = "brown"
cutoff = round(376*0.1)

my_tissue = "Spleen"
my_module = "green"
cutoff = round(101*0.15)
# =======================================================

if(my_tissue == "Liver"){
  my_networkType = "signed hybrid"
  my_power= 5
  my_nrow = 2 # for liver 
  my_height = 5 # for liver 
  
  datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
  net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
} else{my_networkType = "signed hybrid"
my_power= 14
my_nrow = 1 # for spleen
my_height = 2.5 # for spleen 

datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
}

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
  limma_result_of_interest <- limma_list$status} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
  limma_result_of_interest <- limma_list_spleen$cort }


networkdf_STRING <- read_csv(glue("results_networks/cytohubba_STRING_withlogFC_{my_module}_{my_tissue}.csv")) %>% 
  rename(symbol = node_name) %>% 
  select(symbol, MCC) %>% 
  mutate(MCCrank = rank(desc(MCC))) %>% 
  arrange(MCCrank) %>% 
  as.data.frame()

kIN_df <- readRDS(glue("results_RNAseqRDS/kIN_dataframe_{my_tissue}.RDS")) %>%
  filter(module == my_module) %>% 
  arrange(desc(kWithin)) %>% 
  mutate(kINrank = rank(desc(kWithin)))

# ===============================================================
networkdf_STRING %>% 
  head(cutoff) 

kIN_df %>% 
  head(cutoff)


networkdf_STRING %>% 
  head(cutoff) %>%  .$symbol -> x

kIN_df %>% 
  head(cutoff) %>% .$symbol -> y

x[x %in% y] -> overlap
c(x, y) %>% unique() -> both 
overlap
both
`%notin%` <- Negate(`%in%`)
both[both %notin% overlap] -> both_but_overlap

overlap_list[[my_module]] <-overlap
both_list[[my_module]] <-both
bothbutoverlap_list[[my_module]] <-both_but_overlap
# now network plot ===========================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleColors %>% table()
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
modNames = substring(names(MEs), 3)

modules = my_module

# Select module probes
probes = colnames(datExpr)

inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = grcm38$symbol[match(modProbes, grcm38$ensgene)];

selected <- (moduleColors == my_module)
TOM = TOMsimilarityFromExpr(datExpr[,selected],
                            power = my_power,
                            networkType = "signed hybrid")

cyt = exportNetworkToCytoscape(TOM,
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])


set.seed(928)
myweight_cut = median(cyt$edgeData$weight)
cyt$edgeData %>% 
  select(fromNode, toNode, weight) -> edges

cyt$nodeData -> nodes




net <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE) %>% 
  igraph::simplify(., remove.multiple = T) %>% 
  ggnetwork::ggnetwork(.)

`%notin%` <- Negate(`%in%`)
both[both %notin% overlap] -> meh

weight_threshold <- median(net$weight, na.rm = T)

median(net$weight, na.rm = T)

net1 <- net %>% 
  mutate(labelyes = ifelse(name %in% overlap, as.character(name),"")) %>% 
  mutate(labelyes2 = ifelse(name %in% meh, as.character(name),"")) %>% 
  rename(symbol = name) %>%
  mutate(weightx = ifelse(weight < weight_threshold,0,weight)) # change weight cutoff

net2 <- net1 %>% left_join(limma_result_of_interest %>% 
                   select(symbol,logFC) %>% 
                   unique()) %>% 
mutate(log2FC = log2(exp(logFC)))

color_limit = max(abs(net2$log2FC))*(1.05)

mid_rescaler <- function(mid = 0) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}

p <- ggplot(net2,aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey50",alpha=0.2)+
  geom_nodes(aes(color = log2FC,size = weightx),alpha = 0.6,size=8) +
  geom_nodetext(aes( label = labelyes),
                fontface = "bold",
                size =5.5) +
  # geom_nodetext(aes(label = labelyes2), size =4.2) +
  # scale_color_gradient2()+
  scale_color_distiller(palette = "Spectral",limit = c(-color_limit,color_limit),rescaler = mid_rescaler())+
  theme_blank()+
  theme(legend.position = "right")

print(p) 


png(filename = glue("results_figures/hubnetwork_{my_module}_{my_tissue}.png"),
    width = my_width, height = my_height, units = "cm", res = 600)

print(p) 
invisible(dev.off())


# two modules together? =======================================================================
# saveRDS(overlap_list,"results_RNAseqRDS/overlap_list.RDS")
# saveRDS(both_list,"results_RNAseqRDS/both_list.RDS")
# saveRDS(bothbutoverlap_list,"results_RNAseqRDS/bothbutoverlap_list.RDS")

overlap_list <- readRDS("results_RNAseqRDS/overlap_list.RDS")
both_list <- readRDS("results_RNAseqRDS/both_list.RDS")
bothbutoverlap_list <- readRDS("results_RNAseqRDS/bothbutoverlap_list.RDS")

modules = c("pink","cyan")
selected_genes <- c(both_list$pink,both_list$cyan)
overlap <- c(overlap_list$pink, overlap_list$cyan)
bothbutoverlap <- c(bothbutoverlap_list$pink, bothbutoverlap_list$cyan)


selected_genes <- c(both_list$blue)
overlap <- c(overlap_list$blue)
bothbutoverlap <- c(bothbutoverlap_list$blue)


modules = c("brown","green")
selected_genes <- c(both_list$brown,both_list$green)
overlap <- c(overlap_list$brown, overlap_list$green)
bothbutoverlap <- c(bothbutoverlap_list$brown, bothbutoverlap_list$green)

if(my_tissue == "Spleen"){data.frame(symbol = selected_genes) %>% 
    left_join(grcm38 %>% select(symbol,ensgene)) %>% 
    filter(!is.na(ensgene)) -> whenwhen 
} else{
  data.frame(symbol = selected_genes) %>% 
    left_join(grcm38 %>% select(symbol,ensgene)) %>% 
    filter(!is.na(ensgene)) %>% 
    filter(ensgene != "ENSMUSG00000102455") %>% 
    filter(ensgene != "ENSMUSG00000106922") -> whenwhen 
  
}

whenwhen$ensgene -> selected_probes
whenwhen$symbol -> selected_symbols

selected_probes %>% unique() -> mygenename

datexpr_module = datExpr[,mygenename]


modProbes = colnames(datexpr_module)

modGenesymbols = grcm38$symbol[match(modProbes, grcm38$ensgene)]


TOM_hub = TOMsimilarityFromExpr(datexpr_module,
                            power = my_power,
                            networkType = "signed hybrid")


hubprobes = colnames(datexpr_module)
hubGenes = grcm38$symbol[match(hubprobes, grcm38$ensgene)];

dimnames(TOM_hub) = list(hubGenes, hubGenes)
cyt = 
  exportNetworkToCytoscape(TOM_hub,
                           edgeFile = paste("CytoscapeInput-edges-HUB", paste("fourhubs", collapse="-"), ".txt", sep=""),
                           nodeFile = paste("CytoscapeInput-nodes-HUB", paste("fourhubs", collapse="-"), ".txt", sep=""),
                           weighted = TRUE,
                           threshold = 0.02,
                           nodeNames = hubGenes,
                           altNodeNames = hubGenes)





set.seed(928)
myweight_cut = median(cyt$edgeData$weight)
cyt$edgeData %>% 
  select(fromNode, toNode, weight) -> edges

cyt$nodeData -> nodes




net <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE) %>% 
  igraph::simplify(., remove.multiple = T) %>% 
  ggnetwork::ggnetwork(.)


weight_threshold <- median(net$weight, na.rm = T)


net1 <- net %>% 
  mutate(labelyes = ifelse(name %in% overlap, as.character(name),"")) %>% 
  mutate(labelyes2 = ifelse(name %in% bothbutoverlap, as.character(name),"")) %>% 
  rename(symbol = name) %>%
  mutate(weightx = ifelse(weight < weight_threshold,0,weight)) # change weight cutoff

net2 <- net1 %>% left_join(limma_result_of_interest %>% 
                             select(symbol,logFC) %>% 
                             unique()) %>% 
  mutate(log2FC = log2(exp(logFC)))

color_limit = max(abs(net2$log2FC))*(1.05)

mid_rescaler <- function(mid = 0) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}

p <- ggplot(net2,aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(color = "grey50",alpha=0.2)+
  geom_nodes(aes(color = log2FC,size = weightx),alpha = 0.6,size=8) +
  geom_nodetext(aes( label = labelyes),
                fontface = "bold",
                size =5.5) +
  # geom_nodetext(aes(label = labelyes2), size =4.2) +
  # scale_color_gradient2()+
  scale_color_distiller(palette = "Spectral",limit = c(-color_limit,color_limit),rescaler = mid_rescaler())+
  theme_blank()+
  theme(legend.position = "right")



png(filename = glue("results_figures/hubnetwork_{my_module}_{my_tissue}xxxx.png"),
    width = 36, height = 24, units = "cm", res = 300,
    bg = "transparent")

print(p) 
invisible(dev.off())





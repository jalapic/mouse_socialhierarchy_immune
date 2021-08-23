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


# datTraits <- rnaseq_sampleid %>% 
#   filter(sampleID %in% rownames(datExpr)) %>% 
#   left_join(alldata) %>% 
#   mutate(status = ifelse(domgroup == "Alpha", 1,0)) %>% 
#   mutate(pair_status = ifelse(pair_status == "dom", 1,0)) %>% 
#   column_to_rownames("sampleID") %>% 
#   dplyr::select(`B cells_PD09`, `CD4_CD8_ratio_PD09`, `Cytotoxic T_PD09`, `DCs_PD09`,
#     `Helper T_PD09`, `Macro_PD09`, `MLR_PD09`, `Mono_PD09`,
#     `Neutro_PD09`, `NK cells_PD09`, `NLR_PD09`, `T cells_PD09`
#   )


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
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# saveRDS(MEs,glue("results_RNAseqRDS/{my_tissue}_MEs.RDS"))

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
# 
# sizeGrWindow(10,6)
# 
# dev.off() # make sure you do this before AND after
# png(file = glue("results_figures/module_trait_{my_tissue}.png"),
#     width=800, height=890, res = 130)
# 
# # Will display correlations and their p-values
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
#                     signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = glue("{my_tissue}: Module-trait relationships"))
# 
# 
# dev.off()
# # remember - they are correlation, not regression

#=====================================================================================
#
#  Code chunk - Won's code to look at regression 
#
#=====================================================================================



MEs %>%
  rownames_to_column("sampleID") %>%
  left_join(datTraits %>%
              rownames_to_column("sampleID")) %>% 
  mutate(status = ifelse(status ==1, "Alpha", "Subordinate"))-> ME_df





ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(spleen_weight_mg, value, color = status, fill = status))+
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme(legend.position = "top")+
  labs(x = "Spleen weight after 14 days of group housing",
       y = "Module eigengene") -> p1
p1

ggsave(filename = glue("results_figures/{my_tissue}_eigengene_spleenweight_status.png"),
       p1,
       height = my_height, width = 16, dpi = 150)

ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(`Helper T_GD14`, value, color = status, fill = status))+
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme(legend.position = "top")+
  labs(x = "Helper T cell (GD14)",
      y = "Module eigengene") -> p2

p2

ggsave(filename = glue("results_figures/{my_tissue}_eigengene_helperTGD14_status.png"),
       p2,
       height = my_height, width = 16, dpi = 150)



ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  # filter(modules %in% c("MEmagenta","MEcyan","MEgreen")) %>%
  # mutate(status = ifelse(status == 1, "Alpha","Subordinate")) %>%
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme(legend.position = "none") +
  labs(x = "Social status",
       y = "Module eigengene",
       title = glue("{my_tissue}: Module eigengene across social status")) -> p3

p3

ggsave(filename = glue("results_figures/{my_tissue}_eigengene_status_boxplot.png"),
       p3,
       height = my_height, width = 14, dpi = 100)



ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(cort_post, value))+
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme(legend.position = "top")+
  labs(x = "Plasma corticosterone level (GD14)",
       y = "Module eigengene") -> p4

p4

ggsave(filename = glue("results_figures/{my_tissue}_eigengene_cortGD14_status.png"),
       p4,
       height = my_height, width = 16, dpi = 150)



ME_df %>%
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(`Mono_GD14`, value,color = status, fill = status))+
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme(legend.position = "top")+
  labs(x = "Monocyte  (GD14)",
       y = "Module eigengene") -> p5

p5

ggsave(filename = glue("results_figures/{my_tissue}_eigengene_monocyteGD14_status.png"),
       p5,
       height = my_height, width = 16, dpi = 150)



ME_df %>% 
  ggplot(aes(MEpink,MEcyan))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(ME_df$MEcyan,ME_df$MEpink, method = "spearman") # rho -0.86 # p-value = 2.193e-06


ME_df %>% 
  ggplot(aes(MEgreen,MEbrown))+
  geom_point()+
  geom_smooth(method = "lm")

cor.test(ME_df$MEgreen,ME_df$MEbrown, method = "spearman") # rho -0.97 # p-value = 2.168e-06

cor.test(ME_df$MEgreen,ME_df$MEbrown, method = "pearson") # R -0.95 # p-value = 3.708e-12



#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable David's score containing the David's score column of datTrait
status = as.data.frame(datTraits$status);
names(status) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, status, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(status), sep="");
names(GSPvalue) = paste("p.GS.", names(status), sep="");


geneModuleMembership %>% 
  rownames_to_column("ensgene") %>% 
  left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS_all



# to make holistic freaking dataframe ==================================================
# gene_MM_TS_all$module <- moduleColors
# 
# 
# get_mm <- function(x){
#   x$moduleMembership <- x[colnames(x) == paste("MM",x$module,sep = "")] %>% 
#     unlist %>% 
#     as.numeric 
#   xx <- x %>% 
#     dplyr::select(ensgene,module,moduleMembership,GS.trait)
#   return(xx)
# }
# 
# 
# wgcna_whole <- get_mm(gene_MM_TS_all[1,])
# 
# for(i in 2:nrow(gene_MM_TS_all)){
#   wgcna_whole <- rbind(wgcna_whole,get_mm(gene_MM_TS_all[i,]))  
# }
# 
# wgcna_whole %>% 
#   rename(GS.status = GS.trait) -> wgcna_whole
# 
# cort_post = as.data.frame(datTraits$cort_post);
# names(cort_post) = "cort"
# geneTraitSignificance = as.data.frame(cor(datExpr, cort_post, use = "p"));
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
# 
# names(geneTraitSignificance) = paste("GS.", names(cort_post), sep="");
# names(GSPvalue) = paste("p.GS.", names(cort_post), sep="");
# 
# wgcna_whole %>% 
#   left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> wgcna_all
# 
# saveRDS(wgcna_all, glue("results_RNAseqRDS/WGCNA_MM_GS_all_{my_tissue}.RDS"))

# Define cort GD14 (mostly for spleen) ==================================
cort_post = as.data.frame(datTraits$cort_post);
names(cort_post) = "trait"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cort_post, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cort_post), sep="");
names(GSPvalue) = paste("p.GS.", names(cort_post), sep="");



#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
my_trait = "Status"
module = "pink"
module = "blue"
module = "cyan"
module = "grey60"

module_list = c("pink","cyan", "blue", "grey60") # For liver


my_trait = "Corticosterone level (GD14)"
module = "green"
module = "brown"

module_list = c("green", "brown") # For spleen ~ with cort GD14

hub_gene_list = vector('list', length = length(module_list))
names(hub_gene_list) <- module_list


column = match(module, modNames);
moduleGenes = moduleColors==module;
colnames(datExpr)[moduleColors==module] -> module_gene


sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = glue("Gene significance for {my_trait}"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
datExpr[moduleGenes,]

my_MM_threshold = 0.8
my_GS_threshold = 0.2


for (module in module_list){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  colnames(datExpr)[moduleColors==module] -> module_gene
  
  
  grcm38 %>% 
    select(ensgene, symbol, chr, description) %>% 
    filter(ensgene %in% module_gene) -> module_gene_info
  

  
  geneModuleMembership %>% 
    rownames_to_column("ensgene") %>% 
    left_join(geneTraitSignificance %>% rownames_to_column("ensgene")) -> gene_MM_TS
  
  
  
  gene_MM_TS %>% 
    filter(ensgene %in% module_gene) %>% 
    select(ensgene, glue("MM{module}"), GS.trait) %>% 
    filter(abs(GS.trait) >= my_GS_threshold) -> x
  
  x[x[,glue("MM{module}")]>my_MM_threshold,] -> hub_genes
  
  hub_genes %>% 
    left_join(module_gene_info) %>% 
    mutate(moduleName = glue("{module}")) %>% 
    rename(moduleMembership = glue("MM{module}")) -> hub_genes
  hub_gene_list[[module]] <- hub_genes
  
}

hub_gene_list %>% 
  rbindlist() %>% 
  unique() -> hubgenes_df

hubgenes_df %>% 
  filter(symbol!="Sema4d") %>% 
  arrange(desc(GS.trait)) %>% 
  select(-ensgene, -chr) %>% 
  head(15)

write.csv(hubgenes_df, glue("results_tables/hubgene_list_{my_tissue}_cort14.csv"), row.names = F)






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



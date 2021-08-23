S
bio_pcs <- readRDS('results/biological_data_pcs_missingvalue_imputed.RDS')
# principal components of bio data, and behavior data ~ module eigen gene


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


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);





MEs %>%
  rownames_to_column("sampleID") %>%
  left_join(rnaseq_sampleid) %>%
  left_join(bio_pcs) -> ME_pcs_df



ME_pcs_df %>% 
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(Dim.1, value))+
  geom_point()+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme_bw()+
  labs(x = "Principal components 1",
       y = "Module eigengene") -> pp1

pp1
ggsave(filename = glue("results_figures/{my_tissue}_eigengene_bio_PC1.png"),
       pp1,
       height = my_height, width = 14, dpi = 150)


ME_pcs_df %>% 
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(Dim.2, value))+
  geom_point()+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme_bw()+
  labs(x = "Principal components 1",
       y = "Module eigengene") -> pp2

pp2
ggsave(filename = glue("results_figures/{my_tissue}_eigengene_bio_PC2.png"),
       pp2,
       height = my_height, width = 14, dpi = 150)


ME_pcs_df %>% 
  gather(modules,value,2:(1+moduleNumber) ) %>%
  ggplot(aes(Dim.3, value))+
  geom_point()+
  geom_smooth(method = "lm", alpha = 0.2)+
  facet_wrap(~modules, scales = "free_y", nrow = my_nrow)+
  theme_bw()+
  labs(x = "Principal components 1",
       y = "Module eigengene") -> pp3

pp3
ggsave(filename = glue("results_figures/{my_tissue}_eigengene_bio_PC3.png"),
       pp3,
       height = my_height, width = 14, dpi = 150)
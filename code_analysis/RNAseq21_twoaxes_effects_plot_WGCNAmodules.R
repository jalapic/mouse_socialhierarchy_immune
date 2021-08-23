my_tissue = "Liver"
my_module = "pink"

my_tissue = "Liver"
my_module = "cyan"

my_tissue = "Liver"
my_module = "blue"

my_tissue = "Spleen"
my_module = "brown"

my_tissue = "Spleen"
my_module = "green"
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
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
modNames = substring(names(MEs), 3)

modules = my_module

# Select module probes
probes = colnames(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = grcm38$symbol[match(modProbes, grcm38$ensgene)]



limma_list$status %>% 
  filter(symbol %in% modGenes) %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_status = logFC,
         pval_status = P.Value) -> status_df

limma_list$cort %>% 
  filter(symbol %in% modGenes) %>% 
  select(symbol, logFC, P.Value) %>% 
  rename(logFC_cort = logFC,
         pval_cort = P.Value) -> cort_df


df <- full_join(status_df, cort_df) %>% 
  mutate(Significance = ifelse(pval_status<0.05 & pval_cort <0.05, "Both", ifelse(pval_status<0.05 | pval_cort <0.05,"One","None"))) %>% 
  mutate(Significance = factor(Significance, levels =c("Both","One","None")))

max(abs(c(df$logFC_status, df$logFC_cort)))*1.05 -> lim



png(filename = glue("results_figures/twoaxis_{my_tissue}_{my_module}.png"),
    width = 7, height = 8, units = "cm", res = 600)
df %>% 
  ggplot(aes(logFC_status, logFC_cort,color = Significance ))+
  geom_abline(slope = 1, intercept = 0, color ='grey', linetype = 'dashed')+
  geom_hline(yintercept = 0, color ='grey', linetype = 'dashed')+
  geom_vline(xintercept = 0, color ='grey', linetype = 'dashed')+
  geom_point(shape = 21, size = 2)+
  scale_x_continuous(limits = c(-lim,lim))+
  scale_y_continuous(limits = c(-lim,lim))+
  scale_color_brewer(palette = "Dark2")+
  labs(x = "Status effect (Sub <-> Dom)",
       y = "CORT effect (Low <-> High)",
       color = "")+
  theme_bw(base_size = 7)+
  theme(legend.position = "top")
invisible(dev.off())


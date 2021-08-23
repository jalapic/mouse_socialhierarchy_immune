library(CoDiNA)
library(wTO)
h_gene_sets = msigdbr(species = "Mus musculus", category = "H")

my_tissue = "Liver"
my_tissue = "Spleen"

datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  

# subset for the selected hallmark gene set
hallmark = "HALLMARK_INFLAMMATORY_RESPONSE"
my_hallmark = "Inflammatory Response" # for title 
my_pval_cutoff = 0.2
h_gene_sets %>% 
  filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") -> hm_gene_df

hallmark = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
my_hallmark = "Oxidative phosphorylation"
my_pval_cutoff = 0.1
h_gene_sets %>% 
  filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") -> hm_gene_df

hm_gene_df$gene_symbol[hm_gene_df$gene_symbol %in% grcm38$symbol == F]
hm_gene_df$gene_symbol[hm_gene_df$entrez_gene %in% grcm38$entrez == F]



hm_gene_df %>% 
  rename(entrez = entrez_gene) %>%
  # rename(symbol = gene_symbol) %>%
  # left_join(grcm38 %>% 
              # select(symbol, ensgene)) %>% 
  left_join(grcm38 %>% 
              select(entrez, ensgene)) %>% 
  filter(!is.na(ensgene)) %>% 
  .$ensgene %>% unlist -> hm_ensgene

datExpr  %>%  t -> t_datExpr
  
t_datExpr[rownames(t_datExpr) %in% hm_ensgene,] %>% t ->hm_datExpr

str(hm_datExpr)
# now separate for alpha and sub 
rnaseq_sampleid %>% 
  left_join(all_behavior %>% select(subjectID, domgroup)) %>% 
  filter(domgroup == "Alpha" & tissue == my_tissue) %>% 
  .$sampleID -> alphaID

rnaseq_sampleid %>% 
  left_join(all_behavior %>% select(subjectID, domgroup)) %>% 
  filter(domgroup != "Alpha" & tissue == my_tissue) %>% 
  .$sampleID -> subID

hm_datExpr[rownames(hm_datExpr) %in% alphaID,] %>% 
  t %>%  
  as.data.frame %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  select(-ensgene) %>% 
  filter(!is.na(symbol)) %>% 
  filter(!duplicated(symbol)) %>% 
  column_to_rownames('symbol')-> alpha_expression

hm_datExpr[rownames(hm_datExpr) %in% subID,] %>% 
  t %>%  
  as.data.frame %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% select(ensgene, symbol)) %>% 
  select(-ensgene) %>% 
  filter(!is.na(symbol)) %>% 
  filter(!duplicated(symbol)) %>% 
  column_to_rownames('symbol')-> sub_expression

wTO_out_alpha = wTO.fast(Data = alpha_expression, n = 100)
saveRDS(wTO_out_alpha,glue("results_RNAseqRDS/wTO_out_alpha_{my_tissue}_{hallmark}.RDS"))
# wTO_out_alpha <- readRDS(glue("results_RNAseqRDS/wTO_out_alpha_{my_tissue}_{hallmark}.RDS"))

wTO_out_sub = wTO.fast(Data = sub_expression, n = 100)
saveRDS(wTO_out_sub,glue("results_RNAseqRDS/wTO_out_sub_{my_tissue}_{hallmark}.RDS"))
# wTO_out_sub <- readRDS(glue("results_RNAseqRDS/wTO_out_sub_{my_tissue}_{hallmark}.RDS"))

# # explore ============================================================
# wTO_out_alpha %>% 
#   ggplot(aes(-log10(pval)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Alpha wTO unadjusted p-value"))
# 
# wTO_out_sub %>% 
#   ggplot(aes(-log10(pval)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Sub wTO unadjusted p-value"))
# 
# wTO_out_alpha %>% 
#   ggplot(aes(-log10(pval.adj)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Alpha wTO adjusted p-value"))
# 
# wTO_out_sub %>% 
#   ggplot(aes(-log10(pval.adj)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Sub wTO adjusted p-value"))
# 
# 
# # not so pleasant looking  
# 
# wTO_out_alpha %>% 
#   filter(pval.adj <0.05) %>% 
#   nrow()
# 
# wTO_out_sub %>% 
#   filter(pval.adj <0.05) %>% 
#   nrow()
# 
# 
# wTO_out_sub %>% 
#   arrange(pval)
# 
# # ALPHA = subset(wTO_out_alpha, wTO_out_alpha$pval.adj <0.05, select = c('Node.1', 'Node.2', 'wTO'))
# # SUB = subset(wTO_out_sub, wTO_out_sub$pval.adj <0.05, select = c('Node.1', 'Node.2', 'wTO'))
# 
# 
# wTO_out_alpha %>% 
#   arrange(desc(wTO))
# 

# wTO_out_alpha %>% 
#   mutate(pval = replace(pval,pval == 0,0.001)) %>% 
#   mutate(wTO_logpval = sin(wTO)*-log10(pval)) -> wTO_out_alpha
# 
# wTO_out_sub %>% 
#   mutate(pval = replace(pval,pval == 0,0.001)) %>% 
#   mutate(wTO_logpval = sin(wTO)*-log10(pval)) -> wTO_out_sub

# my_cutoff = 0.57 # liver
# my_cutoff = 2.15 # spleen
# 
# wTO_out_alpha %>% 
#   mutate(status = "Alpha") %>% 
#   rbind(wTO_out_sub %>% mutate(status = "Subordinate")) %>% 
#   mutate(pval = replace(pval,pval == 0,0.001)) %>% 
#   mutate(wTO_logpval = sin(wTO)*-log10(pval)) %>% 
#   ggplot(aes(wTO_logpval))+
#   geom_histogram()+ # interesting 
#   geom_vline(xintercept = c(my_cutoff,-(my_cutoff)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: sin(wTO)*-log10(pval)"),
#        x = "sin(wTO)*-log10(pval)")+
#   facet_wrap(~status)

# P.value cutoff 
ALPHA = subset(wTO_out_alpha, wTO_out_alpha$pval < my_pval_cutoff, select = c('Node.1', 'Node.2', 'wTO'))
SUB = subset(wTO_out_sub, wTO_out_sub$pval < my_pval_cutoff, select = c('Node.1', 'Node.2', 'wTO'))


# P.adjust cutoff 
# ALPHA = subset(wTO_out_alpha, wTO_out_alpha$pval.adj < 0.05, select = c('Node.1', 'Node.2', 'wTO'))
# SUB = subset(wTO_out_sub, wTO_out_sub$pval.adj < 0.05, select = c('Node.1', 'Node.2', 'wTO'))


# Run DiffNet
DiffNet = MakeDiffNet(Data = list(ALPHA, SUB), Code = c('ALPHA', 'SUB'))
print(DiffNet) %>% head()

DiffNet %>% 
  filter(Phi == "b")

# Using the median
int_C = quantile(DiffNet$Score_internal, 0.5)
ext_C = quantile(DiffNet$Score_Phi, 0.5)

# quantile
# int_C = quantile(DiffNet$Score_internal, 0.25)
# ext_C = quantile(DiffNet$Score_Phi, 0.75)

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = ext_C, cutoff.internal = int_C)
table(Nodes_Groups$Phi_tilde)

Graph = plot(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, layout = 'layout_components', path = glue('Vis_{my_tissue}_{hallmark}.html'))
# # Using the first and the third quantile
# int_C = quantile(DiffNet$Score_internal, 0.25)
# ext_C = quantile(DiffNet$Score_Phi, 0.75)
# 
# Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = ext_C, cutoff.internal = int_C)
# table(Nodes_Groups$Phi_tilde)



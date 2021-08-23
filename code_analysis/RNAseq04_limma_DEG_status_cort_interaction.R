#dataframes from carpentry/10_get_RNAseq_data.R


my_tissue = "Spleen"
spleen_counts %>% 
  column_to_rownames('ensgene') %>% 
  select(-`G.5`) -> counts

my_tissue = "Liver"
liver_counts %>% 
  column_to_rownames('ensgene') %>% 
  dplyr::select(-`G.5`) -> counts

dim(counts)

# filter_counts = 50
# countData <- counts[rowSums(counts > filter_counts) > round((length(c))*0.9), ]
# dim(countData)

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

d <- d0
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left
  




snames <- colnames(counts) # Sample names
snames


rnaseq_sampleid %>% 
  filter(subjectID!='G.5') %>% 
  filter(tissue == my_tissue) %>% 
  left_join(alldata %>% select(subjectID, domgroup, cort_post)) -> var_info

plotMDS(d, col = ifelse(var_info$domgroup == "Alpha", 10,4))


scale(as.numeric(var_info$cort_post)) -> cort
ifelse(var_info$domgroup == "Alpha", -1,1) -> status

# var_info$domgroup -> status



# ## Model comparison
# # https://rdrr.io/bioc/limma/man/selectmodel.html
# 
# designlist <- list(
#   None=cbind(Int=rep(1,length(cort))),
#   cort=cbind(Int=1,cort=cort),
#   status=cbind(Int=1,status=status), # 0 = alpha, 1 = subordinate
#   Both=cbind(Int=1,cort_status=cort*status),
#   Add=cbind(Int=1,cort=cort,status=status),
#   Full=cbind(Int=1,cort=cort,status=status,cort_status=cort*status)
# )
# out_aic <- selectModel(y, designlist, criterion="aic")
# table(out_aic$pref)
# 
# out_bic <- selectModel(y, designlist, criterion="bic")
# table(out_bic$pref)
# 
# 
# # okkkkkay.. full model then!

mm_full <- model.matrix(~ cort*status)
colnames(mm_full)
head(mm_full)

y <- voom(d, mm_full, plot = T)

# saveRDS(y, glue("results_RNAseqRDS/limma_y_{my_tissue}"))

fit_full <- lmFit(y, mm_full)

#https://mdozmorov.github.io/BIOS567/assets/presentation_diffexpression/DiffExpr_Limma.html


results <- decideTests(fit_full,
                       adjust.method = "none", 
                       p.value = 0.05,
                       lfc = 0)
vennDiagram(results[,-1], include=c("up","down"))

tmp_cort <- contrasts.fit(fit_full, coef = 2) 
tmp_status <- contrasts.fit(fit_full, coef = 3) 
tmp_interaction <- contrasts.fit(fit_full, coef = 4) 

# coef = 1: intercept
# coef = 2: cort
# coef = 3: status
# coef = 4: cort:status (interaction)

tmp_cort <- eBayes(tmp_cort) 
tmp_status <- eBayes(tmp_status)
tmp_interaction <- eBayes(tmp_interaction)


topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val,description)  -> limma_cort

topTable(tmp_status, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val,description) -> limma_status

topTable(tmp_interaction, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val,description) -> limma_interaction


# limma_list <- list()
# limma_list$cort <- limma_cort
# limma_list$status <- limma_status
# limma_list$interaction <- limma_interaction
# 
# 
# saveRDS(limma_list,glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
# 


limma_interaction %>% 
  filter(symbol == "Nr3c1")

limma_cort %>% 
  filter(symbol == "Nr3c1")


limma_status %>% 
  filter(symbol == "Fkbp5")

limma_cort %>% 
  filter(symbol == "Fkbp5")

limma_cort %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(abs(logFC))) -> limma_cort_DEG

limma_status %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(abs(logFC))) -> limma_status_DEG

limma_interaction %>% 
  filter(P.Value <0.05) %>% 
  arrange(desc(abs(logFC))) -> limma_interaction_DEG

nrow(limma_interaction_DEG)

head(limma_cort_DEG,20)
head(limma_status_DEG,20)
head(limma_interaction_DEG,20)

my_symbol = "Mup20"
my_symbol = "Mup3"
my_symbol = "Mup17"
my_symbol = "Mup16"
my_symbol = "Mup14"
my_symbol = "Gcgr"
my_symbol = "Igf1"
my_symbol = "Cdk9"
my_symbol = "Gpd2"
my_symbol = "Mup21"
my_symbol = "Mup2"
my_symbol = "Mup5"
my_symbol = "Mup18"
my_symbol = "Ly6g6c"
my_symbol = "Terb1"
my_symbol = "Hrg"
my_symbol = "Mmd2" # monocyte to macrophage differentiation-associated 2
my_symbol = "Nr3c1"
my_symbol = "Fkbp5"

grcm38 %>% 
  filter(symbol == my_symbol) %>% 
  .$ensgene %>% 
  unique() -> goi


y$E[goi,] %>%
  as.data.frame(col.names = c("Expression")) %>% 
  rename(Expression = ".") %>%  
  rownames_to_column("subjectID") %>% 
  left_join(var_info) %>% 
  ggplot(aes(cort_post,Expression, color = domgroup, fill = domgroup)) + 
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  labs(title = glue("{my_tissue}: {my_symbol}"),
      x = "Corticosterone ng/ml (GD14)",
      y = "Gene expression")+
  theme(legend.position = "top")






# ================================================
my_symbols = c("Mup20","Mup3", "Mup17","Mup16", "Mup14", "Gcgr", "Igf1", "Cdk9", "Gpd2", "Mup21", "Mup2", "Mup5", "Mup18", "Ly6g6c", "Terb1", "Mmd2")


grcm38 %>% 
  filter(symbol == my_symbols) %>% 
  .$ensgene %>% 
  unique() -> goi


y$E[goi,] %>%
  as.data.frame(col.names = c("Expression")) %>% 
  t()
  rename(Expression = ".") %>%  
  rownames_to_column("subjectID") %>% 
  left_join(var_info) %>% 
  ggplot(aes(cort_post,Expression, color = domgroup, fill = domgroup)) + 
  geom_point(size = 3, shape = 21, alpha = 0.3)+
  geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  # labs(title = glue("{my_tissue}: {my_symbol}"),
  #      x = "Corticosterone ng/ml (GD14)",
  #      y = "Gene expression")+
  theme(legend.position = "top")+
  facet_wrap( ~ symbol)


# ========================================================
write.csv(limma_cort_DEG,glue("results_tables/{my_tissue}_limma_cort_DEG.csv"),row.names = F)
write.csv(limma_status_DEG,glue("results_tables/{my_tissue}_limma_status_DEG.csv"), row.names = F)
write.csv(limma_interaction_DEG,glue("results_tables/{my_tissue}_limma_interaction_DEG.csv"),row.names = F)

## CORT only =====================================
# mm <- model.matrix(~ cort)
# head(mm)
# y <- voom(d, mm, plot = T)
# 
# fit <- lmFit(y, mm)
# tmp <- contrasts.fit(fit, coef = 2) # test "pH" coefficient
# tmp <- eBayes(tmp)
# 
# top.table <- topTable(tmp, sort.by = "P", n = Inf)
# 
# top.table %>% 
#   rownames_to_column('ensgene') %>% 
#   left_join(grcm38) %>%
#   select(symbol,logFC,P.Value,adj.P.Val,t,description) -> limma_result
# 
# head(limma_result,20)
# 
# limma_result %>% 
#   filter(P.Value <0.05) -> limma_DEG
# 
# limma_result %>% 
#   filter(P.Value <0.05) %>% 
#   filter(logFC >0) -> limma_DEG_upregulated_with_cort
# 
# 
# limma_result %>% 
#   filter(P.Value <0.05) %>% 
#   filter(logFC <0) -> limma_DEG_downregulated_with_cort
# 
# nrow(limma_DEG)
# nrow(limma_DEG_upregulated_with_cort)
# nrow(limma_DEG_downregulated_with_cort)





q.dl.liver <- readRDS("results_RNAseqRDS/limma_eFDR_Liver.RDS")
q.dl.spleen <- readRDS("results_RNAseqRDS/limma_eFDR_Spleen.RDS")


head(q.dl.liver)
head(q.dl.spleen)















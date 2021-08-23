



# my data 

my_tissue = "Spleen"
spleen_counts %>% 
  column_to_rownames('ensgene') %>% 
  select(-`G.5`) -> dlNorm


my_tissue = "Liver"
liver_counts %>%
  column_to_rownames('ensgene') %>%
  dplyr::select(-`G.5`) -> dlNorm


rnaseq_sampleid %>% 
  filter(subjectID!='G.5') %>% 
  filter(tissue == my_tissue) %>% 
  left_join(alldata %>% select(subjectID, domgroup, cort_post)) -> var_info

colnames(dlNorm)
# var_info$domgroup %>% 
#   factor(.,levels = c("Alpha","Subordinate")) -> group.dl 
ifelse(var_info$domgroup == "Alpha", 1,-1) -> group.dl
var_info$cort_post %>% 
  scale %>% 
  as.numeric -> cort.dl 

d = apply(dlNorm, 2, as.numeric)
d0= DGEList(d, group = group.dl)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)


design.dl <- model.matrix(~ group.dl*cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl, glue("results_RNAseqRDS/limma_vdl_{my_tissue}"))
# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(i in 1 : R){
  print(paste("Starting on Permutation", i))

  # Randomize the traits
  cort.dl.rand = sample(cort.dl)
  group.dl.rand = sample(group.dl)

  # Model
  design.dl.rand = model.matrix(~group.dl.rand*cort.dl.rand)
  colnames(design.dl.rand) <- mycolnames

  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
}

q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

for(i in 1 : R){
  print(paste("Calculating Permutation", i))

  temp = p.dl.rand[[i]]

  for(c in 1 : 4){
    for(r in 1 : nrow(p.dl.limma)){
      if(temp[r, c] <= p.dl.limma[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)

saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_2000_tworand.RDS"))



q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_2000_tworand.RDS"))


head(round(p.dl.limma,5))
head(q.dl)
hist(q.dl-p.dl.limma)

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()



# ========================================================================
# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))

tmp_status <- contrasts.fit(efit.dl, coef = 2) # group
tmp_cort <- contrasts.fit(efit.dl, coef = 3) # cort
tmp_interaction <- contrasts.fit(efit.dl, coef = 4) # interaction

limma_list <- list()

topTable(tmp_status, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$status


topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cort


topTable(tmp_interaction, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$interaction


saveRDS(limma_list,glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))



# example from Joyce and Becca (Hoffman Lab) ==============================
# # Expression values
# dlNorm = read.csv("dlNorm.csv", row.names = 1)
# 
# # Group
# group.dl = c("hour", "hour", "day", "week", "day", "hour", "hour", "week",
#              "week", "week", "week", "day", "day", "hour", "hour", "day",
#              "day", "day", "day", "hour")
# group.dl = factor(group.dl, levels = c("hour", "day", "week"))
# 
# # Create a DGEList
# dge.dl = apply(dlNorm, 2, as.numeric)
# dge.dl = DGEList(dge.dl, group = group.dl)
# 
# # p-values
# p.dl.limma = read.csv("p.dl.limma.csv")
# 
# # Traits (aggression scores)
# aggr.dl = c(11.00000, 5.00000, 31.00000, 6.00000, 15.00000, 0.00000, 25.00000,
#             30.00000, 26.14286, 32.00000, 30.00000, 10.00000, 16.00000,
#             13.00000, 14.00000, 19.00000, 17.60000, 0.00000, 21.00000, 2.00000)
# 
# # How many random sampling
# R = 5000
# 
# p.dl.rand = list()
# 
# for(i in 1 : R){
#   print(paste("Starting on Permutation", i))
#   
#   # Randomize the traits
#   aggr.dl.rand = sample(aggr.dl)
#   
#   # Model
#   design.dl.rand = model.matrix(~group.dl + aggr.dl.rand)
#   colnames(design.dl.rand) = c("hour", "day", "week", "aggr")
#   
#   # Calculate p-values based on randomized traits
#   v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
#   vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
#   efit.dl.rand = eBayes(vfit.dl.rand)
#   p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
# }
# 
# q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))
# 
# for(i in 1 : R){
#   print(paste("Calculating Permutation", i))
#   
#   temp = p.dl.rand[[i]]
#   
#   for(c in 1 : 4){
#     for(r in 1 : nrow(p.dl.limma)){
#       if(temp[r, c] <= p.dl.limma[r, c]){
#         q.dl[r, c] = q.dl[r, c] + 1
#       }
#     }
#   }
# }
# 
# q.dl = q.dl / R
# colnames(q.dl) = c("hour", "day", "week", "aggr")
# q.dl = as.data.frame(q.dl)
# 



# my data 

my_tissue = "Spleen"
spleen_counts %>% 
  column_to_rownames('ensgene') %>% 
  select(-`G.5`) -> dlNorm


my_tissue = "Liver"
liver_counts %>%
  column_to_rownames('ensgene') %>%
  dplyr::select(-`G.5`) -> dlNorm
# How many random sampling
R = 5000
# ==============================================================================
limma_list <- list()

rnaseq_sampleid %>% 
  filter(subjectID!='G.5') %>% 
  filter(tissue == my_tissue) %>% 
  left_join(alldata %>% select(subjectID, domgroup, cort_post)) -> var_info

colnames(dlNorm)
var_info$domgroup %>%
  factor(.,levels = c("Alpha","Subordinate")) -> group.dl
# ifelse(var_info$domgroup == "Alpha", 1,-1) -> group.dl
var_info$cort_post %>% 
  scale %>% 
  as.numeric -> cort.dl 

d = apply(dlNorm, 2, as.numeric)
d0= DGEList(d)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff) 
dge.dl <- d0[-drop,]
dim(dge.dl)


# status ==============================================================================
design.dl <- model.matrix(~ group.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]


# p.dl.rand = vector('list',length = R)
# 
# for(i in 1 : R){
#   print(paste("Starting on Permutation", i))
# 
#   # Randomize the traits
#   group.dl.rand = sample(group.dl)
# 
#   # Model
#   design.dl.rand = model.matrix(~group.dl.rand)
#   colnames(design.dl.rand) <- mycolnames
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
#   for(c in 1 : 2){
#     for(r in 1 : nrow(p.dl.limma)){
#       if(temp[r, c] <= p.dl.limma[r, c]){
#         q.dl[r, c] = q.dl[r, c] + 1
#       }
#     }
#   }
# }
# 
# q.dl = q.dl / R
# colnames(q.dl) <- mycolnames
# q.dl = as.data.frame(q.dl)
# row.names(q.dl) <- rownames(dge.dl)
# 
# saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_status.RDS"))

q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_status.RDS"))
hist(p.dl.limma - q.dl)
head(q.dl)

# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))

tmp_status <- eBayes(efit.dl)

topTable(efit.dl, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) %>% 
  arrange(desc(abs(logFC))) %>% 
  mutate(logFC = (-1)*logFC) -> limma_list$status

# 
# limma_list_liver$status %>% 
#   arrange(desc(abs(logFC))) %>%  
#   filter(symbol == "Myc")


# cort ========================================================================

design.dl <- model.matrix(~ cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]


# p.dl.rand = vector('list',length = R)
# 
# for(i in 1 : R){
#   print(paste("Starting on Permutation", i))
# 
#   # Randomize the traits
#   cort.dl.rand = sample(cort.dl)
# 
#   # Model
#   design.dl.rand = model.matrix(~cort.dl.rand)
#   colnames(design.dl.rand) <- mycolnames
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
#   for(c in 1 : 2){
#     for(r in 1 : nrow(p.dl.limma)){
#       if(temp[r, c] <= p.dl.limma[r, c]){
#         q.dl[r, c] = q.dl[r, c] + 1
#       }
#     }
#   }
# }
# 
# q.dl = q.dl / R
# colnames(q.dl) <- mycolnames
# q.dl = as.data.frame(q.dl)
# row.names(q.dl) <- rownames(dge.dl)
# 
# saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_cort.RDS"))

q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_cort.RDS"))

# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))


topTable(efit.dl, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort

hist(q.dl-p.dl.limma)


limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}_second.RDS"))













# NO interACTION ==============================================================================
design.dl <- model.matrix(~ group.dl+cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = T)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]


p.dl.rand = vector('list',length = R)

for(i in 1 : R){
  print(paste("Starting on Permutation", i))
  
  # Randomize the traits 
  group.dl.rand = sample(group.dl)
  cort.dl.rand = sample(cort.dl)
  
  # Model
  design.dl.rand = model.matrix(~group.dl.rand + cort.dl.rand)
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
  
  for(c in 1 : 3){
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

saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_noint_tworand.RDS"))

q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_tissue}_cutoff5_{R}_noint_tworand.RDS"))
hist(p.dl.limma - q.dl)
hist(p.dl.limma[2])
hist(q.dl[[2]])
# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- rownames(efit.dl$t)
nrow(q.dl)
nrow(efit.dl$t)
sum(duplicated(row.names(efit.dl$coefficients)))



tmp_status <- contrasts.fit(efit.dl, coef = 2) # group
tmp_cort <- contrasts.fit(efit.dl, coef = 3) # cort


topTable(tmp_status, sort.by = "P", n = Inf, adjust ="fdr") %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) %>% 
  arrange(desc(abs(logFC))) %>% 
  mutate(logFC = (-1)*logFC) -> limma_list$status

topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) %>% 
  arrange(desc(abs(logFC))) -> limma_list$cort

saveRDS(limma_list, glue("results_RNAseqRDS/limma_{my_tissue}_add.RDS"))

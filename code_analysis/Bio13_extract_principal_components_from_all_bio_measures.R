# PCA - with complete cases, too few samples 

pca_df <- alldata %>% 
  mutate(status = ifelse(domgroup == "Alpha", 1,0)) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", 1,0)) %>% 
  column_to_rownames('subjectID') %>% 
  select(
    GD01_BW, GD14_bw, spleen_weight_mg, AGD, 
    cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
    `B cells_GD14`, CD4_CD8_ratio_GD14, `Cytotoxic T_GD14`, `DCs_GD14`,
    `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
    `Neutro_GD14`, `NK cells_GD14`, `T cells_GD14`
  )


prcomp(pca_df %>% na.omit,
       scale = TRUE,
       center = TRUE) -> res.pca

res.pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var ='subjectID') %>% 
  left_join(alldata %>% select(subjectID, domgroup)) %>% 
  ggplot(aes(PC1, PC2,color = domgroup))+
  geom_point(size = 3, alpha = 0.5) + 
  scale_color_manual(values = c("purple4", "green","orange")) 



# Using special packages to impute the missing values ===================

# install.packages("VIM")
# install.packages('missMDA')
# install.packages("FactoMineR")
library(VIM)
library(missMDA)
library(FactoMineR)

res<-summary(aggr(pca_df, sortVar=TRUE))$combinations


nb <- estim_ncpPCA(pca_df,method.cv = "Kfold", verbose = FALSE) # estimate the number of components from incomplete data
#(available methods include GCV to approximate CV)
nb$ncp #3

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

res.comp <- imputePCA(pca_df, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

imp <- cbind.data.frame(res.comp$completeObs,cormat_df$status)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = nb$ncp, graph=FALSE)

plot(res.pca, hab=12, lab="quali")


res.pca$ind$coord %>% 
  as.data.frame() %>% 
  rownames_to_column('subjectID') %>% 
  left_join(alldata %>% select(subjectID,domgroup)) %>% 
  ggplot(.,aes(`Dim.1`, `Dim.2`,color = domgroup))+
  geom_point(size = 3, alpha = 0.5) +
  theme(legend.position = "top")+
  labs(x = "PC1", y = "PC2",
       color = "Social status")+
  ggtitle("PCA plots with all biological data (except RNAseq data), missing data imputed") -> bio_pc_plot

bio_pc_plot

ggsave("results_figures/biological_data_PCA.png", 
       bio_pc_plot,
       width = 7,
       height = 5.5, 
       dpi = 150
       )


res.pca$ind$coord %>%
  as.data.frame() %>%
  rownames_to_column('subjectID') %>%
  left_join(alldata %>% select(subjectID,domgroup)) %>%
  saveRDS('results_statRDS/biological_data_pcs_missingvalue_imputed.RDS')

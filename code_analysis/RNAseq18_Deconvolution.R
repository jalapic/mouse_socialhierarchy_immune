# https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html

# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")

library(immunedeconv)

set_cibersort_binary("C:/Users/PSYC-wl8856/Desktop/Github/CIBERSORT.R")
set_cibersort_mat("C:/Users/PSYC-wl8856/Desktop/Github/LM22.txt")


res = deconvolute(immunedeconv::dataset_racle$expr_mat, "quantiseq")
knitr::kable(res, digits=2)


res_epic = deconvolute(immunedeconv::dataset_racle$expr_mat, "epic")
knitr::kable(res_epic, digits=2)

# CIBERSORT - special case 

res_cibersort = deconvolute(immunedeconv::dataset_racle$expr_mat, "cibersort")
knitr::kable(res_cibersort, digits=2)



# how is the input data looking like?
immunedeconv::dataset_racle$expr_mat -> x

glimpse(x)





# Dataset conversion
## rownames are expected to be in HGNC gene symbols 
## OR, The Bioconductor ExpressionSet works 
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertMouseGeneList <- function(x){

  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  
  if(length(x)==length(humanx)){return(humanx)} # this if else was the charm! 
  else{return(NA)}
}



# mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
#       values = rownames(spleen_exp), 
#       mart = mart)     


# my current RNAseq data, let's do spleen first 
spleen_exp <- readRDS("results_RNAseqRDS/datExpr_Spleen.RDS") %>% 
  t()

head(spleen_exp)
rownames(spleen_exp)

spleen_exp %>% 
  as.data.frame() %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>% dplyr::select(ensgene, symbol)) %>% 
  .$symbol -> symbol_to_convert

convertMouseGeneList('Cdc45')
convertMouseGeneList('Narf')

# spleen_exp %>% 
#   as.data.frame() %>% 
#   rownames_to_column('ensgene') %>% 
#   left_join(grcm38 %>% dplyr::select(ensgene, symbol)) %>% 
#   mutate(hgnc = convertMouseGeneList(symbol)) -> converted

spleen_exp_hgnc <- vector('character',length = length(symbol_to_convert))
length(symbol_to_convert)
for (i in 1:99){
  spleen_exp_hgnc[i] <- convertMouseGeneList(symbol_to_convert[i])
}

saveRDS(spleen_exp_hgnc, "spleen_exp_hgnc.RDS" )



res = deconvolute(spleen_exp_hgnc, "quantiseq")
knitr::kable(res, digits=2)



saveRDS(x,"spleen_deconv_input.RData")

write.csv(x,"spleen_deconv_input.csv")

head(x)

any(is.na(x))
tail(x)
sum(is.na(row.names(x)))

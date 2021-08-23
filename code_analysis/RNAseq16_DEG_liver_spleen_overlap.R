limma_list_spleen <- readRDS("results_RNAseqRDS/limma_Spleen.RDS") %>% 
  map(~select(.,-description))
limma_list_liver <- readRDS("results_RNAseqRDS/limma_Liver.RDS") %>% 
  map(~select(.,-description))

# Status 

limma_list_spleen$status %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> x1


limma_list_liver$status %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> y1

inner_join(x1,y1) %>% 
  unique() -> status_both


status_both %>% 
  left_join(wgcna_all_liver) %>% 
  arrange(module)

# Cort
limma_list_spleen$cort %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> x2


limma_list_liver$cort %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> y2

inner_join(x2,y2) %>% 
  unique()



# Interaction
limma_list_spleen$interaction %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> x3


limma_list_liver$interaction %>% 
  filter(P.Value <0.05) %>% 
  filter(abs(logFC) < log2(1.2)) %>% 
  select(symbol) -> y3

inner_join(x3,y3) %>% 
  unique()

# Protein-protein network with STRING database 

# http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_BS32010/Workshops/BS32010Workshop4.html
# http://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf 




string_db <- STRINGdb$new( version="11", species=10090, # for mouse
              score_threshold=200, input_directory="")
# as of 11/23/2020, version 11.0 is available on STRING

limma_list_spleen <- readRDS("results_RNAseqRDS/limma_Spleen.RDS") %>% 
  map(~select(.,-description))
limma_list_liver <- readRDS("results_RNAseqRDS/limma_Liver.RDS") %>% 
  map(~select(.,-description))

wgcna_all_spleen <- readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_Spleen.RDS") %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) %>% select(-ensgene)
wgcna_all_liver <- readRDS("results_RNAseqRDS/WGCNA_MM_GS_all_Liver.RDS") %>% 
  left_join(grcm38 %>% select(ensgene,symbol)) %>% select(-ensgene)

module_list = c("pink","cyan", "blue", "grey60") # For liver
module_list = c("green", "brown") # For spleen ~ with cort GD14


# Liver 
my_tissue = "Liver"

my_module = "pink"
my_GS_threshold = 0.0

my_module = "cyan"
my_GS_threshold = 0.0



wgcna_all_liver %>% 
  filter(module == my_module) %>% 
  filter(abs(GS.status) >my_GS_threshold) %>% 
  .$symbol -> ppi_symbol
  
limma_list_liver$status %>% 
  filter(symbol %in% ppi_symbol) %>% 
  arrange(P.Value) %>%  
  left_join(grcm38 %>% dplyr::select(symbol,entrez)) %>% 
  rename(pvalue = P.Value,
         gene = entrez) %>% 
  dplyr::select(pvalue, logFC, gene) -> ppi_gene


limma_list_liver$status %>% 
  filter(symbol %in% ppi_symbol) %>% 
  select(symbol,logFC) %>% 
write.table(.,
            glue("results_networks/for_ppi_{my_module}_{my_tissue}.txt"),
            row.names = F, quote = F, sep = "\t", col.names = F)





# Spleen 
my_tissue = "Spleen"

my_module = "brown"
my_MM_threshold = 0.0


my_module = "green"
my_MM_threshold = 0.0


wgcna_all_spleen %>% 
  filter(module == my_module) %>% 
  filter(moduleMembership >my_MM_threshold) %>% 
  .$symbol -> ppi_symbol

limma_list_spleen$cort %>% # because spleen modules are more associated with cort 
  filter(symbol %in% ppi_symbol) %>% 
  arrange(P.Value) %>%  
  left_join(grcm38 %>% dplyr::select(symbol,entrez)) %>% 
  rename(pvalue = P.Value,
         gene = entrez) %>% 
  dplyr::select(pvalue, logFC, gene) -> ppi_gene

limma_list_spleen$cort %>% 
  filter(symbol %in% ppi_symbol) %>% 
  select(symbol,logFC) %>% 
  write.table(.,
              glue("results_networks/for_ppi_{my_module}_{my_tissue}.txt"),
              row.names = F, quote = F, sep = "\t", col.names = F)





# https://rdrr.io/bioc/STRINGdb/src/R/rstring.R
# map it ====================================================

my_mapped <- string_db$map(ppi_gene, "gene", removeUnmappedRows = TRUE )

hits <- my_mapped$STRING_id
# firsthit <- my_mapped$STRING_id[1]
# string_db$plot_network( hits )

# filter by p-value and add a color column
# (i.e. green down-regulated gened and red for up-regulated genes)
my_mapped_pval05 <- string_db$add_diff_exp_color(subset(my_mapped, pvalue<0.05),
                                                       logFcColStr="logFC" )
# post payload information to the STRING server
payload_id <- string_db$post_payload(my_mapped_pval05$STRING_id,
                                     colors=my_mapped_pval05$color )
# display a STRING network png with the "halo"

png(filename = glue("results_figures/PPI_{my_module}module_{my_tissue}.png"),
    width = 2400, height = 2400, res =600)

string_db$plot_network( hits, payload_id=payload_id )

dev.off()

annotations <- string_db$get_annotations( hits )
head(annotations, n=20)
string_db$get_graph()



wgcna_all_liver %>% 
  filter(module == my_module) %>% 
  arrange(desc(moduleMembership)) 

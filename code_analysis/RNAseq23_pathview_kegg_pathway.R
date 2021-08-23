# BiocManager::install('pathview')
library(pathview)

my_tissue = "Spleen"
my_tissue = "Liver"
my_logFC_threshold = 0.0


limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS")) %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

limma_status_DEG <- limma_list$status
limma_cort_DEG <- limma_list$cort
limma_interaction_DEG <- limma_list$interaction

getKEGG <- function(limma_df,my_showCategory = 10){
  kegg_df <- limma_df %>%
    left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
    mutate(entrez = as.character(entrez)) %>% 
    dplyr::select(entrez, logFC) %>% 
    filter(!is.na(entrez)) %>% 
    arrange(desc(logFC))
  
  kegg <- enrichKEGG(gene    = kegg_df$entrez,
                     organism     = 'mmu',
                     pvalueCutoff = 0.05)
  
  fortify(
    kegg,
    showCategory = my_showCategory,
    by = "GeneRatio",
    split = NULL,
    includeAll = TRUE ) %>% 
    arrange(desc(GeneRatio))  -> temp1
 
  return(temp1)
  
}



my_kegg_list <- limma_list %>% map(~getKEGG(.,my_showCategory = 200)) %>% 
  map(~select(., ID, Description))

my_kegg_list

# Spleen ==========================================================
# status
pathview_df <- limma_status_DEG %>% 
  left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(!is.na(entrez)) %>% 
  filter(!duplicated(entrez)) 
logFC <- pathview_df$logFC
names(logFC) <- pathview_df$entrez

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu00190",  # Oxidative phosphorylation
                   species = "mmu")

# cort 
pathview_df <- limma_cort_DEG %>% 
  left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(!is.na(entrez)) %>% 
  filter(!duplicated(entrez)) 
logFC <- pathview_df$logFC
names(logFC) <- pathview_df$entrez

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04062",  # Chemokine signaling pathway
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04662",  # B cell receptor signaling pathway
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04660",  # T cell receptor signaling pathway
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04659",  # Th17 cell differentiation
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04658", #    Th1 and Th2 cell differentiation
                   species = "mmu")

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04064",  # NF-kappa B signaling pathway
                   species = "mmu")



# Liver ==========================================================
# status
pathview_df <- limma_status_DEG %>% 
  left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(!is.na(entrez)) %>% 
  filter(!duplicated(entrez)) 
logFC <- pathview_df$logFC
names(logFC) <- pathview_df$entrez

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu00061",  # FATTY ACID BIOSYNTHESIS
                   species = "mmu")
pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04979", # CHOLESTEROL METABOLISM
                   species = "mmu")
pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu00140", # STERIOD HORMONE BIOSYNTHESIS
                   species = "mmu")
pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04935", #  Growth hormone synthesis, secretion and action
                   species = "mmu")

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu00190",  # Oxidative phosphorylation
                   species = "mmu")

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04935",  # growth hormone 
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu01212",  # fatty acid metabolism 
                   species = "mmu")



pv.out <- pathview(gene.data = logFC, 
                  pathway.id = "mmu04660",  # T cell receptor signaling
                  species = "mmu")

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04658",  # Th1 and Th2 cell differentiation
                   species = "mmu")


pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04064",  # NF-kappa B signaling pathway
                   species = "mmu")

 

# cort 
pathview_df <- limma_cort_DEG %>% 
  left_join(grcm38 %>% dplyr::select(symbol, entrez)) %>% 
  filter(!is.na(entrez)) %>% 
  filter(!duplicated(entrez)) 
logFC <- pathview_df$logFC
names(logFC) <- pathview_df$entrez

pv.out <- pathview(gene.data = logFC, 
                   pathway.id = "mmu04612",  # Antigen processing and presentation
                   species = "mmu")

v.out <- pathview(gene.data = logFC, 
                  pathway.id = "mmu04660",  # T cell receptor signaling
                  species = "mmu")










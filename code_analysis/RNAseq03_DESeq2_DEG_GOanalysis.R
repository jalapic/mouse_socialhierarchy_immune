data(geneList, package="DOSE")

my_tissue = "Liver"

my_tissue = "Spleen"

go_df <- read_csv(glue('results_tables/{my_tissue}_DEG_table.csv')) %>% 
  left_join(grcm38 %>% dplyr::select(ensgene, entrez)) %>% 
  # filter(log2FoldChange > 0) %>%
  # filter(log2FoldChange < 0) %>%
  mutate(entrez = as.character(entrez)) %>% 
  dplyr::select(tissue,entrez) %>% 
  filter(!is.na(entrez)) 



ggo <- enrichGO(gene = go_df$entrez ,
                OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                readable = T,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.10)


dev.off()
png(glue("results_figures/enrichGO_{my_tissue}.png"), width = 8, height = 4, units = 'in', res = 300)
clusterProfiler::dotplot(ggo, orderBy = "Count")+
  ggtitle(glue("{my_tissue}"))+
  theme(plot.title = element_text(size = 12))

dev.off()


clusterProfiler::goplot(ggo)


ggo %>% 
  filter(Description == "response to wounding") %>%
  as.data.frame() %>% 
  rename(symbols = geneID) %>% 
  .$symbols %>% 
  (str_split(.,"/") %>% 
  as.data.frame(col.names = "symbol") %>% )
  left_join(grcm38 %>% select(symbol,description))

ggo %>% 
  filter(Description == "wound healing") %>% 
  as.data.frame() %>% 
  rename(symbols = geneID) %>% 
  .$symbols %>% 
  str_split(.,"/") %>% 
  as.data.frame(col.names = "symbol") %>% 
  left_join(grcm38 %>% select(symbol,description))


ggo %>% 
  filter(Description == "steroid metabolic process") %>%
  as.data.frame() %>% 
  rename(symbols = geneID) %>% 
  .$symbols %>% 
  str_split(.,"/") %>% 
  as.data.frame(col.names = "symbol") %>% 
  left_join(grcm38 %>% select(symbol,description))


ggo %>% 
  filter(Description == "fatty acid metabolic process") %>%
  as.data.frame() %>% 
  rename(symbols = geneID) %>% 
  .$symbols %>% 
  str_split(.,"/") %>% 
  as.data.frame(col.names = "symbol") %>% 
  left_join(grcm38 %>% select(symbol,description))







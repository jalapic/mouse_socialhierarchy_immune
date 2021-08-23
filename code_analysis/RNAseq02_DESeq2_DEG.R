my_tissue = 'Liver'
my_tissue = 'Spleen'




# adjusted pval and fold change threshold to 'define' DEG
adj_pval_cutoff = 0.1
threshold = 0.9  # 90% of the subject
LFC_threshold = 0.2
pvalue_threshold = 0.05 # didn't use it for now


dds <- readRDS(glue("RNAseq_RDS/dds_{my_tissue}.RDS"))

temp <- results(dds, contrast=c("domgroup","Alpha","Subordinate")) %>% 
  as.data.frame() %>% 
  rownames_to_column('ensgene') %>% 
  as_tibble() %>% 
  select(ensgene, log2FoldChange, pvalue, padj) %>% 
  mutate(contrast = "Alpha - Subordinate") %>% 
  filter(abs(log2FoldChange)>=LFC_threshold) %>%
  filter(pvalue <=pvalue_threshold) %>%
  mutate(tissue = my_tissue) %>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  left_join(grcm38 %>%  
              mutate(description = gsub("\\[.*?\\]",'',description)) %>% 
              select(ensgene, symbol, chr, biotype, description)) %>% 
  unique()

# write.csv(temp, glue('tables/{my_tissue}_DEG_table.csv'),row.names = F)

temp %>% 
  head(25) %>% 
  as.data.frame() %>% 
  select(symbol,log2FoldChange,pvalue,padj, description) %>% 
  arrange(desc(abs(log2FoldChange)))

temp %>% 
  filter(grepl('Mup',symbol)) %>% 
  as.data.frame() %>% 
  select(symbol,log2FoldChange,pvalue,padj, description) %>% 
  arrange(as.numeric(gsub('Mup','',symbol)))


temp %>% 
  filter(grepl('inter',description)) %>% 
  as.data.frame() %>% 
  select(symbol,log2FoldChange,pvalue,padj, description) %>% 
  arrange(abs(log2FoldChange))




# MUP plot - liver===============================================================
results(dds, contrast=c("domgroup","Alpha","Subordinate")) %>% 
  as.data.frame() %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38 %>%  
              mutate(description = gsub("\\[.*?\\]",'',description)) %>% 
              select(ensgene, symbol, chr, biotype, description)) %>% 
  filter(grepl('Mup',symbol)) %>% 
  arrange(as.numeric(gsub('Mup','',symbol))) -> mup_gene_df 

mup_gene_df %>% 
  .$ensgene %>% 
  unique() -> mup_gene_list

mup_gene_df %>% 
  .$symbol %>% 
  unique() -> mup_symbol_list



grcm38 %>% 
  filter(symbol %in% mup_gene_list) %>% 
  select(ensgene,symbol) %>% 
  arrange(as.numeric(gsub('Mup','',symbol))) %>% 
  unique()-> goi
goi %>% as.data.frame()


goi$ensgene -> goi_ensgene
goi_ensgene %in% names(dds)

goi %>% filter(ensgene %in% names(dds))

coldata

tcounts <- t(log2((counts(dds[mup_gene_list, ], replaced=FALSE)+.5))) %>%
  as.data.frame( ) %>% 
  rownames_to_column( "sampleID") %>% 
  gather(ensgene, expression, (ncol(.)-length(mup_gene_list)+1):ncol(.)) %>% 
  left_join(grcm38 %>% select(symbol, ensgene)) %>% 
  left_join(coldata) %>% 
  mutate(symbol = factor(symbol, levels = mup_symbol_list))



ggplot(tcounts, aes(symbol, expression,color=domgroup,fill= domgroup)) + 
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, 
                 alpha = 0.3,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.6)) +
  labs(x="", 
       y="Expression (log normalized counts)",
       fill = "Status",
       color= "Status")+
  scale_color_manual(values = c('purple4','orange'))+
  scale_fill_manual(values = c('purple4','orange'))+
  annotate('rect', xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 7.5, xmax = 8.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 13.5, xmax = 14.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 15.5, xmax = 16.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 17.5, xmax = 18.5, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  annotate('rect', xmin = 19.5, xmax = Inf, ymin = -Inf, ymax = Inf,fill = "grey", alpha = 0.2) +
  # scale_color_viridis(discrete = T,option = "C")+
  # scale_fill_viridis(discrete = T,option = "C")+
  theme_bw()+
  theme(legend.position = "top",
        panel.grid.major = element_blank()) -> p


print(p)


# ggsave(p, filename = "figures/Mups_liver.png",height = 5,width = 10,dpi = 200)


# tick plot =====================================================================
# just curious the pattern of DEG in liver ~ spleen, alpha vs. subordinate 











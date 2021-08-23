get_DEG_results <- function(dds, LFC_threshold = 0.2, pvalue_threshold = 0.05) # I won't set padj yet 
{
  
  temp <- results(dds, contrast=c("domgroup","Alpha","Subordinate")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "Alpha - Subordinate") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>% 
    filter(pvalue <=pvalue_threshold) 
  
  
  temp %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}

# functions used in RNAseq14_extract_gene_info_DEG_WGCNA.R
name_any_gene <- function(my_gene)
{
  print(grcm38 %>% 
          filter(symbol == my_gene) %>% 
          .$description)
  print("Liver")
  
  print("~ Status")
  print(limma_list_liver$status %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print("~ CORT (GD14)")
  print(limma_list_liver$cort %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print("~ CORT:Status")
  print(limma_list_liver$interaction %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print(wgcna_all_liver %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol))
  
  print("Spleen")
  
  print("~ Status")
  print(limma_list_spleen$status %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print("~ CORT (GD14)")
  print(limma_list_spleen$cort %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print("~ CORT:Status")
  print(limma_list_spleen$interaction %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol, -description))
  
  print(wgcna_all_spleen %>% 
          filter(symbol ==my_gene) %>% 
          select(-symbol))
  
}

plot_any_gene <- function (my_gene){
  print(grcm38 %>% 
          filter(symbol == my_gene) %>% 
          .$description)
  
  grcm38 %>% 
    filter(symbol == my_gene) %>% 
    .$ensgene %>% 
    unique() -> goi
  
  if (goi %in% rownames(y_liver$E) == F)
  {temp_p1 <- ggplot()+ggtitle("No detectable gene expression in liver")}
  else {
    y_liver$E[goi,] %>%
      as.data.frame(col.names = c("Expression")) %>% 
      rename(Expression = ".") %>%  
      rownames_to_column("subjectID") %>% 
      left_join(alldata %>% 
                  select(subjectID, cort_post,domgroup)) %>% 
      ggplot(aes(cort_post,Expression, color = domgroup, fill = domgroup)) + 
      geom_point(size = 3, shape = 21, alpha = 0.3)+
      geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
      scale_color_manual(values = c("purple4","orange"))+
      scale_fill_manual(values = c("purple4","orange"))+
      labs(title = glue("Liver: {my_gene}"),
           x = "Corticosterone ng/ml (GD14)",
           y = "Gene expression")+
      theme_bw()+
      theme(legend.position = "top") -> temp_p1}
  
  if ( goi %in% rownames(y_spleen$E) == F)
  {temp_p2 <- ggplot()+ggtitle("No detectable gene expression in spleen")}
  else {
    y_spleen$E[goi,] %>%
      as.data.frame(col.names = c("Expression")) %>% 
      rename(Expression = ".") %>%  
      rownames_to_column("subjectID") %>% 
      left_join(alldata %>% 
                  select(subjectID, cort_post,domgroup)) %>% 
      ggplot(aes(cort_post,Expression, color = domgroup, fill = domgroup)) + 
      geom_point(size = 3, shape = 21, alpha = 0.3)+
      geom_smooth(method = "lm", alpha = 0.2, size = 1.2)+
      scale_color_manual(values = c("purple4","orange"))+
      scale_fill_manual(values = c("purple4","orange"))+
      labs(title = glue("Spleen: {my_gene}"),
           x = "Corticosterone ng/ml (GD14)",
           y = "Gene expression")+
      theme_bw()+
      theme(legend.position = "top") -> temp_p2}
  grid.arrange(temp_p1,temp_p2)}


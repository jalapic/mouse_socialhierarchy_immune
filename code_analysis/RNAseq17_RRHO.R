dds_list <- list()

dds_list$Spleen <- readRDS(glue("RNAseq_RDS/dds_spleen.RDS"))
dds_list$Liver <- readRDS(glue("RNAseq_RDS/dds_liver.RDS"))


dds_list %>% 
  map(~get_DEG_results(., LFC_threshold = 0, pvalue_threshold = 1)) %>%  # so pretty much getting everything
  map2_df(.,names(.), ~mutate(.x,tissue = .y)) %>% 
  left_join(grcm38 %>%  select(ensgene, symbol, chr, biotype, description)) %>% 
  unique() -> rrho_df

rrho_df %>% 
  split(.$tissue) %>% 
  map(~select(.,ensgene)) %>% 
  map(unique) -> temp_list

temp_list[[1]] %>% 
  inner_join(temp_list[[2]]) %>% 
  unlist -> genes_expressed_in_both_tissue

length(genes_expressed_in_both_tissue) 

rrho_df %>% 
  filter(ensgene %in% genes_expressed_in_both_tissue) -> rrho_df #limit to genes that are expressed in both liver and spleen


set.seed(1521)

# Regular RRHO

RRHO_tissue_1 = unique(rrho_df$tissue)[1]
RRHO_tissue_2 = unique(rrho_df$tissue)[2]

my_contrast = unique(rrho_df$contrast)

rrho_df %>% 
  filter(contrast == my_contrast) %>% 
  filter(tissue == RRHO_tissue_1) %>% 
  mutate(DDE = -log(pvalue)*log(2^(log2FoldChange))) %>%  
  rename(geneID = ensgene) %>% 
  select(geneID, DDE) %>% 
  filter(!is.na(geneID)) %>% 
  as.data.frame()-> gene_list1

rrho_df %>% 
  filter(contrast == my_contrast) %>% 
  filter(tissue == RRHO_tissue_2) %>% 
  mutate(DDE = -log(pvalue)*log(2^(log2FoldChange))) %>%  
  rename(geneID = ensgene) %>% 
  select(geneID, DDE) %>% 
  filter(!is.na(geneID)) %>% 
  as.data.frame() -> gene_list2



# Regular RRHO 
RRHO_obj <-  RRHO::RRHO(gene_list1, gene_list2, 
                        alternative = "enrichment", 
                        BY = T,
                        log10.ind = T)

# dev.off()
# 
# png(glue('results_figures/RRHO_{RRHO_tissue_1}_{RRHO_tissue_2}.png'),
#     width = 3.5, height = 3.5, units = "in", res = 300)
# try(RRHO2_heatmap(RRHO_obj))
# mtext(glue("{RRHO_tissue_2}"), 2, 0.5)
# mtext(glue("{RRHO_tissue_1}"), 1, 0.5)
# 
# dev.off()



# RRHO2 package : http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html

RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, 
                               labels = c(RRHO_tissue_1, RRHO_tissue_2), 
                               method = "hyper",
                               log10.ind=T)

# dev.off()
# 
# png(glue('results_figures/RRHO2_{RRHO_tissue_1}_{RRHO_tissue_2}.png'),
#     width = 3.5, height = 3.5, units = "in", res = 300)
# RRHO2_heatmap(RRHO_obj)
# mtext(glue("{RRHO_tissue_2}"), 2, 0.5)
# mtext(glue("{RRHO_tissue_1}"), 1, 0.5)
# 
# dev.off()




# RRHO package plot is not working for me. Gotta make my own...?

# ggplot2 
png(glue('results_figures/RRHO_{RRHO_tissue_1}_{RRHO_tissue_2}_ggplot.png'),
    width = 5, height = 3.8, units = "in", res = 100)

melted_hypermat <- reshape2::melt(RRHO_obj$hypermat)

ggplot(data = melted_hypermat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_viridis(na.value = 'white',option = "C")+
  # theme_void()+ 
  labs(x = RRHO_tissue_1,
       y = RRHO_tissue_2,
       fill = "-log10(P-value)")+
  theme(
    line =               element_blank(),
    rect =               element_blank(),
    text =               element_text(
      family = "", face = "plain",
      colour = "black", size = 11,
      lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
      margin = margin(), debug = FALSE
    ),
    axis.text =          element_blank(),
    axis.title.x =         element_text(size = rel(2),vjust = 3),
    axis.title.y =         element_text(size = rel(2),vjust = -0.5,hjust = 0.5,angle = 90),
    axis.ticks.length =  unit(0, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    legend.box =         NULL,
    legend.key.size =    unit(1.2, "lines"),
    legend.position =    "right",
    legend.text =        element_text(size = rel(0.8)),
    legend.title =       element_text(hjust = 0),
    panel.ontop =        FALSE,
    panel.spacing =      unit(11/2, "pt"),
    plot.margin =        unit(c(0, 0, 0, 0), "lines"),
    plot.title =         element_text(
      size = rel(1.2),
      hjust = 0, vjust = 1,
      margin = margin(t = 11/2)
    ),
    plot.title.position = "panel",
    complete = TRUE
  )

dev.off()



# or just work around the function 
maximum = NULL 
minimum = NULL 
colorGradient = NULL


hypermat <- RRHO_obj$hypermat
labels <- c(RRHO_tissue_1, RRHO_tissue_2)
method <- RRHO_obj$method

#
hypermat %>% str
58^2
hypermat[hypermat %>% is.na]


minimum <- min(hypermat, na.rm = TRUE)
maximum <- max(hypermat, na.rm = TRUE)
jet.colors <- colorRampPalette(c("#00007F", "blue", 
                                 "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                 "#FF7F00", "red", "#7F0000"))
colorGradient <- jet.colors(101)

color.bar <- function(lut, min, max = -min, nticks = 11, 
                      ticks = seq(min, max, len = nticks), title = ""){
  scale <- (length(lut) - 1)/(max - min)
  plot(c(0, 10), c(min, max), type = "n", bty = "n", 
       xaxt = "n", xlab = "", yaxt = "n", 
       ylab = "")
  mtext(title, 2, 2.3, cex = 0.8)
  axis(2, round(ticks, 0), las = 1, cex.lab = 0.8)
  for (i in 1:(length(lut) - 1)) {
    y <- (i - 1)/scale + min
    rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
  }}

image(hypermat,  
      col = colorGradient, 
      axes = FALSE)
mtext(labels[2], 2, 0.5)
mtext(labels[1], 1, 0.5)

atitle <- ifelse(RRHO_obj$log10.ind, "-log10(P-value)", 
                 "-log(P-value)")
color.bar(colorGradient, min = min(0, minimum), max = maximum, 
          nticks = 6, title = atitle)



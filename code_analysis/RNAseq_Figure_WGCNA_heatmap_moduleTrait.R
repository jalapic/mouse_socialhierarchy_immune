my_tissue = "Liver"
my_width = 11.6
my_height = 12.5
my_width2 = 17
my_height2 = 13.68


my_tissue = "Spleen"
my_width = 7.8
my_height = 12.5
my_width2 = 17
my_height2 = 6.8


modnum <- readRDS(glue("results_RNAseqRDS/modnum_{my_tissue}.RDS"))

if(my_tissue == "Liver"){
  my_networkType = "signed hybrid"
  my_power= 5
  net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
} else{my_networkType = "signed hybrid"
my_power= 14
net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
}


moduleColors = labels2colors(net$colors)
moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")
modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))

lm_result_all <- readRDS("results_RNAseqRDS/lm_result_list.RDS") %>% 
  rbindlist()

lm_result_all %>% 
  mutate(tissue = ifelse(grepl("Liver",module_name),"Liver","Spleen")) %>% 
  filter(tissue == my_tissue) %>% 
  mutate(module = gsub(".*_", "", module_name)) %>% 
  mutate(point_shape = ifelse(Estimate > 0, "P","N")) -> df

df$module <- factor(df$module,levels=rev(unique(modnum$module[order(modnum$count)])))

lm_result_all %>% 
  select(trait) %>% 
  unique() %>% 
  mutate(timepoint = ifelse(grepl("14",trait),"GD14",
                            ifelse(grepl("post",trait),"GD14",
                                   ifelse(grepl("spleen_weight",trait),"GD14",
                                          ifelse(trait == "status","GD14","PD09"))))) %>% 
  arrange(trait) -> info_temp

new_label = c("AGD",
              "B cells","B cells",
              "CD4/CD8","CD4/CD8",
              "Corticosterone","Corticosterone",
              "Fkbp5 CpG1","Fkbp5 CpG2","Fkbp5 CpG3",
              "Cytotoxic T","Cytotoxic T",
              "DCs","DCs",
              "GD01 BW","GD14 BW",
              "Helper T", "Helper T", 
              "Lymphoid DCs","Lymphoid DCs",
              "Macrophages","Macrophages",
              "MLR","MLR",
              "Monocytes","Monocytes", 
              "Myeloid DCs","Myeloid DCs",
              "Neutrophils","Neutrophils",
              "NK cells","NK cells",
              "NLR","NLR",
              "Pair status","Spleen weight","Status",
              "T cells","T cells")

info_temp %>% cbind(new_label) -> info


info$new_label = factor(info$new_label, 
                        ordered = T,
                        levels = rev(c(
                          "AGD","Pair status",
                          "Status",
                          "GD01 BW","GD14 BW",
                          "Spleen weight","Corticosterone",
                          "Fkbp5 CpG1","Fkbp5 CpG2","Fkbp5 CpG3",
                          "Macrophages","Monocytes","Neutrophils",
                          "NK cells",
                          "Lymphoid DCs","Myeloid DCs", "DCs",
                          "B cells","T cells",
                          "Helper T","Cytotoxic T",
                          "CD4/CD8","MLR","NLR")))


c("Status",
  "Corticosterone",
  "GD14 BW",
  "Spleen weight",
  "Fkbp5 CpG1","Fkbp5 CpG2","Fkbp5 CpG3",
  "Macrophages","Monocytes","Neutrophils",
  "NK cells",  "Lymphoid DCs","Myeloid DCs",
  "B cells","T cells",
  "Helper T","Cytotoxic T","MLR","NLR") -> GD14_level

df %>% 
  left_join(info) -> dfx


dfxx <- dfx %>% 
  filter(timepoint== "GD14") %>% 
  mutate(new_label = factor(new_label, levels = rev(GD14_level))) %>% 
  mutate(my_alpha = ifelse(`Pr(>|t|)` < 0.05, 1, 0)) %>% 
  mutate(my_alpha2 = 1-`Pr(>|t|)` ) %>% 
  filter(module != "grey") %>% 
  filter(new_label != "CD4/CD8") %>% 
  filter(new_label != "DCs")

# for liver 

# saveRDS(dfxx, "results_RNAseqRDS/module_heatmap_df_Liver.RDS")

# if liver,===================================================
# ggplot(data = dfxx %>%  filter(my_alpha2>=0.95),
#        aes(x = module, y =new_label, fill = Estimate))+
#   geom_tile(color = "white",size =1 )+
#   geom_point(data = NULL, aes( x= "pink",y = c("Fkbp5 CpG3")),pch = NA)+
#   scale_y_discrete(limits = rev(GD14_level))+
#   scale_fill_distiller(palette = "PiYG", limit = c(-my_limit,my_limit))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(legend.position = "right", 
#         legend.title = element_text(size =rel(0.75)),
#         legend.text = element_text(size =rel(0.7)),
#         legend.key.size = unit(0.3, 'cm'),
#         axis.ticks = element_blank(),
#         panel.border = element_blank(),
#         # panel.background = element_blank(),
#         axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1)),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2,size = rel(1.5))) -> traitmodule
# 
# 
# 
# png(filename = glue("results_figures/traitmodule_{my_tissue}_newcol2.png"),
#     width = my_width, height = my_height, units = "cm", res = 600)
# 
# traitmodule
# 
# invisible(dev.off())


# ROTATE ==================================================================

dfxx %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(dfxx$module)) -> mehh
mehh[mehh != "grey"] -> y_limit

vlines = c(1:18)+0.495
hlines = c(1:16)+0.51
ggplot(data = dfxx %>%  filter(my_alpha2>=0.95),
       aes(y = module, x =new_label, fill = Estimate))+
  geom_tile(color = "white",size =1 )+
  geom_vline(xintercept = vlines, color = "grey",size = 0.2)+
  geom_hline(yintercept = hlines, color = "grey",size = 0.2)+
  # geom_point(data = NULL, aes( y= "pink",x = c("Fkbp5 CpG3")),pch = NA)+
  scale_x_discrete(limits = GD14_level, position = 'top')+
  scale_y_discrete(limits = y_limit)+
  scale_fill_distiller(palette = "PiYG",
                       limit = c(-my_limit,my_limit),
                       breaks = c(-0.1,0,0.1))+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=0.5,vjust=0.2,size = rel(1)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> traitmodule2


traitmodule2


# 
# png(filename = glue("manuscript/figures/fig4/traitmodule_Liver.png"),
#     width = 16.1, height = 14.5, units = "cm", res = 600)
# 
# print(traitmodule2)
# 
# invisible(dev.off())

# module color and size 

modnum %>% 
  filter(module != "grey") %>% 
  rename(my_alpha2 = count) %>% 
  mutate(module = factor(module, levels = y_limit)) %>% 
  mutate(new_label = "Module size    ") -> modnumx

ggplot(modnumx, aes(y = module,x = new_label, fill = module))+
  geom_tile(color = "white",size =1 )+
  scale_fill_manual(values = y_limit)+
  scale_x_discrete( position = 'top')+
  scale_y_discrete(limits = y_limit, position = "right")+
  geom_text(data = modnumx,  aes(label = my_alpha2))+
  geom_text(data = modnumx %>% filter(module %in% c("black","blue","midnightblue")), 
            aes(label = my_alpha2),
            color = "white")+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> modnumplot
# 
# png(filename = "manuscript/figures/fig4/modulecolorsize_Liver.png",
#     width = 5.8, height = 14.5, units = "cm", res = 600)
# print(modnumplot)
# dev.off()

# 
# cowplot::plot_grid(modnumplot+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#                    traitmodule2+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
#                    rel_widths = c(0.2,0.8),
#                    lign = "h", ncol = 2)
# 



#if spleen, renaming spleen WGCNA module color names ===============================

showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}

nb.cols <- 40
standardColors(n = nb.cols)
showCols(standardColors(n = nb.cols), bg = "white")

newcol <- data.frame(module = factor(c("turquoise","blue","grey","brown","yellow","green","red")),
                     newspleencolor = factor(c("   royalblue    ","steelblue", "grey", "orange",
                                               "     violet     ",
                                               "darkgreen", "darkred")))


dfxxx <- dfxx %>% left_join(newcol) %>% 
  mutate(newspleencolor = factor(newspleencolor,
                                 levels =c("   royalblue    ","steelblue", "grey", "orange",
                                           "     violet     ",
                                           "darkgreen", "darkred")))



# saveRDS("results_RNAseqRDS/module_heatmap_df_Spleen.RDS")
# PLOT ======================================

# dfxxx %>% 
#   filter(my_alpha > 0.95) %>% 
#   .$Estimate %>% 
#   abs() %>% 
#   max() -> my_limit
#   
# 
# ggplot(data = dfxxx %>%  filter(my_alpha2>=0.95),
#        aes(x = newspleencolor, y =new_label, fill = Estimate))+
#   geom_tile(color = "white",size =1 )+
#   geom_point(data = NULL, aes( x= "pink",y = c("Fkbp5 CpG3")),pch = NA)+
#   scale_y_discrete(limits = rev(GD14_level))+
#   scale_fill_distiller(palette = "PiYG", limit = c(-my_limit,my_limit))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(legend.position = "right", 
#         legend.title = element_text(size =rel(0.75)),
#         legend.text = element_text(size =rel(0.7)),
#         legend.key.size = unit(0.3, 'cm'),
#         axis.ticks = element_blank(),
#         panel.border = element_blank(),
#         # panel.background = element_blank(),
#         axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1)),
#         axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2,size = rel(1))) -> traitmodule
# 
# 
# 
# 
# 
# png(filename = glue("results_figures/traitmodule_{my_tissue}_newcol2.png"),
#     width = my_width, height = my_height, units = "cm", res = 600)
# 
# traitmodule
# 
# invisible(dev.off())


# ROTATE ==================================================================

dfxxx %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(dfxxx$newspleencolor)) -> mehh
mehh[mehh != "grey"] -> y_limit

vlines = c(1:18)+0.495
hlines = c(1:5)+0.51
ggplot(data = dfxxx %>%  filter(my_alpha2>=0.95),
       aes(y = newspleencolor, x =new_label, fill = Estimate))+
  geom_tile(color = "white",size =1 )+
  geom_vline(xintercept = vlines, color = "grey",size = 0.2)+
  geom_hline(yintercept = hlines, color = "grey",size = 0.2)+
  # geom_point(data = NULL, aes( y= "pink",x = c("Fkbp5 CpG3")),pch = NA)+
  scale_x_discrete(limits = GD14_level, position = 'top')+
  scale_y_discrete(limits = y_limit)+
  scale_fill_distiller(palette = "PiYG",
                       limit = c(-my_limit,my_limit),
                       breaks = c(-0.1,0,0.1))+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=0.5,vjust=0.2,size = rel(1)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> traitmodule2_spleen



traitmodule2_spleen


# 
# png(filename = glue("manuscript/figures/fig4/traitmodule_Spleen.png"),
#     width = 16.1, height = 8, units = "cm", res = 600)
# 
# traitmodule2
# 
# invisible(dev.off())
# 

 # module color and size 


modnum %>% 
  filter(module != "grey") %>% 
  left_join(newcol) %>% 
  mutate(module = newspleencolor) %>% 
  rename(my_alpha2 = count) %>% 
  mutate(module = factor(module, levels = y_limit)) %>% 
  mutate(new_label = "Module size    ") -> modnumx

ggplot(modnumx, aes(y = module,x = new_label, fill = module))+
  geom_tile(color = "white",size =1 )+
  scale_fill_manual(values = y_limit)+
  scale_x_discrete( position = 'top')+
  scale_y_discrete(limits = y_limit, position = "right")+
  geom_text(data = modnumx,  aes(label = my_alpha2))+
  geom_text(data = modnumx %>% filter(module %in% c("darkred","darkgreen")), 
            aes(label = my_alpha2),
            color = "white")+
  labs(x = "", y = "")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none", 
        legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1))) -> modnumplot_spleen

modnumplot_spleen
# 
# png(filename = "manuscript/figures/fig4/modulecolorsize_Spleen.png",
#     width = 5, height = 8, units = "cm", res = 600)
# print(modnumplot_spleen)
# dev.off()
# 

png(filename = "manuscript/figures/fig4/both.png",
    width = 20, height = 21.7, units = "cm", res = 600)

cowplot::plot_grid(
  modnumplot+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  traitmodule2+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
  modnumplot_spleen+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  traitmodule2_spleen+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  rel_widths = c(0.2,0.8),
  rel_heights = c(0.668,0.332),
  align = "h", ncol = 2)

dev.off()



png(filename = "manuscript/figures/fig4/both2.png",
    width = 19.3, height = 23, units = "cm", res = 600)

cowplot::plot_grid(
  modnumplot+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  traitmodule2+theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm")), 
  modnumplot_spleen+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  traitmodule2_spleen+theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm")),
  rel_widths = c(0.12,0.88),
  rel_heights = c(0.67,0.33),
  align = "h", ncol = 2)

dev.off()

#===============================================













png(filename = glue("results_figures/traitmodule_legend_newcol.png"),
    width = 16.1, height = 14.5, units = "cm", res = 600)


get_legend(traitmodule2) -> leg
grid.arrange(leg)


invisible(dev.off())


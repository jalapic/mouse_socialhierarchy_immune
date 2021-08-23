my_tissue = "Liver"
my_tissue = "Spleen"



c("Module size",
  "Status",
  "Corticosterone",
  "GD14 BW",
  "Spleen weight",
  "Fkbp5 CpG1","Fkbp5 CpG2","Fkbp5 CpG3",
  "Macrophages","Monocytes","Neutrophils",
  "NK cells",  "Lymphoid DCs","Myeloid DCs",
  "B cells","T cells",
  "Helper T","Cytotoxic T","MLR","NLR") -> GD14_level




modnum <- readRDS(glue("results_RNAseqRDS/modnum_{my_tissue}.RDS"))
heatmap_df <- readRDS(glue("results_RNAseqRDS/module_heatmap_df_{my_tissue}.RDS"))
# Liver ===========================================================

heatmap_df %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit

rev(levels(heatmap_df$module)) -> mehh
mehh[mehh != "grey"] -> y_limit


modnum %>% 
  filter(module != "grey") %>% 
  rename(my_alpha2 = count) %>% 
  mutate(module = factor(module, levels = y_limit)) %>% 
  mutate(new_label = "Module size") -> modnumx

ggplot(modnumx, aes(y = module,x = new_label, fill = module))+
  geom_tile(color = "white",size =1 )+
  scale_fill_manual(values = y_limit)+
  scale_x_discrete( position = 'top')+
  scale_y_discrete(limits = y_limit)+
  geom_text(data = modnumx,  aes(label = my_alpha2))+
  geom_text(data = modnumx %>% filter(module %in% c("black","blue","midnightblue")), aes(label = my_alpha2), color = "white")+
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
        axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1.5))) -> modnumplot

png(filename = "results_figures/modulesize2_Liver.png",
    width = 6, height = 14.5, units = "cm", res = 600)
print(modnumplot)
dev.off()


  
heatmap_df %>% 
  full_join(modnum) -> heatmap_dfx


heatmap_dfx %>% 
  filter(new_label == "Module size")
glimpse(heatmap_dfx)
head(heatmap_dfx)
View(heatmap_dfx)

vlines = c(1:18)+0.495
hlines = c(1:16)+0.51
ggplot(data = heatmap_df %>%  filter(my_alpha2>=0.95& my_alpha <2),
       aes(y = module, x =new_label, fill = Estimate))+
  geom_tile(color = "white",size =1 )+
  # geom_text(data = modnum, aes(module,new_label,label = my_alpha2))+
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
  theme(legend.title = element_text(size =rel(1.2)),
        legend.text = element_text(size =rel(1.)),
        legend.position = "none",
        legend.key.size = unit(.9, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1.5))) -> traitmodule2


traitmodule2



png(filename = "results_figures/traitmodule_Liver_rotated_newcol.png",
    width = 17, height = 14.5, units = "cm", res = 600)

traitmodule2

invisible(dev.off())


# Spleen ===========================================================

heatmap_df %>% 
  filter(my_alpha > 0.95) %>% 
  .$Estimate %>% 
  abs() %>% 
  max() -> my_limit
rev(levels(heatmap_df$newspleencolor)) -> mehh
mehh[mehh != "grey"] -> y_limit

vlines = c(1:18)+0.495
hlines = c(1:5)+0.51
ggplot(data = heatmap_df %>%  filter(my_alpha2>=0.95),
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
  # theme_bw()+
  theme(legend.position = "right", 
        legend.title = element_text(size =rel(0.75)),
        legend.text = element_text(size =rel(0.7)),
        legend.key.size = unit(0.3, 'cm'),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "grey"),
        # panel.background = element_blank(),
        # panel.grid.minor = element_line(color = 'grey50'),
        panel.grid = element_blank(),
        axis.text.y = element_text(hjust=0.95,vjust=0.2,size = rel(1)),
        axis.text.x = element_text(angle = 45,hjust=0.08,vjust=0.15,size = rel(1.5))) -> traitmodule2


traitmodule2



png(filename = glue("results_figures/traitmodule_Liver_rotated_newcol.png"),
    width = 17, height = 7.3, units = "cm", res = 600)

traitmodule2

invisible(dev.off())


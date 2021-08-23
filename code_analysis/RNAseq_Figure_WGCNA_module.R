my_networkType = "signed hybrid"
my_tissue = "Liver"
my_power= 5
my_nrow = 2 # for liver 
my_height = 4.5
my_mod_height = 6 # for liver 


my_networkType = "signed hybrid"
my_tissue = "Spleen"
my_power= 14
my_nrow = 1 # for spleen
my_height = 2.1
my_mod_height = 3 # for spleen 


datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
net <- readRDS(glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
dim(datExpr)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

moduleNumber = length(unique(moduleColors))

rownames(datExpr)
colnames(alldata)
alldata$subjectID
datTraits <- rnaseq_sampleid %>% 
  mutate(subjectID = glue("{LETTERS[as.numeric(str_sub(tissueID,1,3))-104]}.{str_sub(tissueID,5)}")) %>% 
  filter(sampleID %in% rownames(datExpr)) %>% 
  left_join(alldata) %>% 
  mutate(status = ifelse(domgroup == "Alpha", 1,0)) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", 1,0)) %>% 
  column_to_rownames("sampleID") %>% 
  dplyr::select(
    GD01_BW, GD14_bw, spleen_weight_mg, AGD, pair_status,
    status,ds, glicko_rank, win_prop, loss_prop, ind_cert,
    preem_flee_perc, flee_perc, sub_pos_perc, freeze_perc, 
    despotism, cort_pre, cort_post, CpG_1_GD14, CpG_2_GD14, CpG_3_GD14,
    `B cells_GD14`, `CD4_CD8_ratio_GD14`, `Cytotoxic T_GD14`, `DCs_GD14`,
    `Helper T_GD14`, `Macro_GD14`, `MLR_GD14`, `Mono_GD14`,
    `Neutro_GD14`, `NK cells_GD14`, `NLR_GD14`, `T cells_GD14`
  )

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

MEs %>%
  rownames_to_column("sampleID") %>%
  left_join(datTraits %>%
              rownames_to_column("sampleID")) %>% 
  mutate(status = ifelse(status ==1, "Alpha", "Subordinate"))-> ME_df



# Module gene numbers =============================================

moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq) -> modnum 
colnames(modnum) <- c("module","count")
str(modnum)

saveRDS(modnum, glue("results_RNAseqRDS/modnum_{my_tissue}.RDS"))


# Plot - for liver =================================================


modnum$module <- factor(modnum$module,levels=unique(modnum$module[order(modnum$count)]))
mycol <- modnum %>% 
  filter(module != "grey") %>% 
  .$module %>% as.character() # for some reason this has to be character again to be ordered properly in the figure...!! 


modnum %>%
  filter(module != "grey") %>% 
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.2, size = 1.4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycol)+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "", 
       title = glue("{my_tissue}: Module size")) -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_tissue}.png"),
    width = 6, height = my_height, units = "cm", res = 600)
temp_p
invisible(dev.off())


modnum %>%
  filter(module != "grey") %>% 
  mutate(module = factor(module, levels = rev(levels(module)))) %>% 
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = 1.1, size = 1.4, angle = 180)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = rev(mycol))+
  theme_minimal(base_size = 5.5)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 180, margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "") -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_tissue}_upsidedown_newcol2.png"),
    width = 6, height = my_height, units = "cm", res = 600)
temp_p
invisible(dev.off())


#  Module gene numbers for spleen - renaming color names =====================================
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
                     newspleencolor = factor(c("royalblue","steelblue", "grey", "orange","violet","darkgreen", "darkred")))

modnumx <- modnum %>% 
  left_join(newcol) %>% 
  filter(module != "grey") %>% 
  mutate(newspleencolor = factor(newspleencolor,
                                 levels = rev(c("royalblue","steelblue", "orange","violet","darkgreen", "darkred"))))
  
mycolx <- modnumx$newspleencolor %>% as.character()
modnumx %>%
  ggplot(aes(y = newspleencolor, x = count, fill = newspleencolor))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = -0.2, size = 1.4)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = mycolx)+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "", 
       title = glue("{my_tissue}: Module size")) -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_tissue}_newcol.png"),
    width = 6, height = my_height, units = "cm", res = 600)
temp_p
invisible(dev.off())


modnumxx <- modnumx %>% 
  mutate(newspleencolor = factor(newspleencolor,levels = c("royalblue","steelblue", "orange","violet","darkgreen", "darkred")))
modnumxx %>%
  ggplot(aes(y = newspleencolor, x = count, fill = newspleencolor))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = 1.1, size = 1.4, angle = 180)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = rev(mycolx))+
  theme_minimal(base_size = 6)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 180, margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "") -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_tissue}_upsidedown_newcol.png"),
    width = 6, height = my_height, units = "cm", res = 600)
temp_p
invisible(dev.off())




# Liver - pink and cyan ME ~ status====================================
standard_size = 5


colnames(ME_df)
ME_df %>% 
  ggplot(aes(status,MEpink, color = status, fill = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 jitter.color = NA,
                  alpha = 0.4,
                 jitter.size = 2.2,
                   jitter.height = 0.02, jitter.width = 0.1, 
                 errorbar.draw = TRUE,
                   position = position_dodge(0.85)) +
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 10)+
  theme(legend.position = "none",
        # axis.title.x = element_blank(),
        plot.title = element_text(size = 10)) +
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Liver-Pink module") -> temp_p2

png(filename = glue("results_figures/MEpink_{my_tissue}_biggerfont.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
temp_p2
invisible(dev.off())

ME_df %>% 
  ggplot(aes(status,MEcyan, color = status, fill = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 jitter.color = NA,
                 alpha = 0.4,
                 jitter.size = 2.2,
                 jitter.height = 0.02, jitter.width = 0.1, 
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 10)+
  theme(legend.position = "none",
        # axis.title.x = element_blank(),
        plot.title = element_text(size = 10)) +
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Liver-Cyan module") -> temp_p3

  
png(filename = glue("results_figures/MEcyan_{my_tissue}_biggerfont.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
temp_p3
invisible(dev.off())

# Liver - blue  ME ~ CORT =================================

ME_df %>% 
  ggplot(aes(cort_post,MEblue))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(size = 1.6, shape = 21, alpha = 0.6,size = 2.2)+
  theme_bw(base_size = 10)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10)) +
  labs(x = "Plasma CORT (GD14)",
       y = "Module eigengene",
       title = "Liver-Blue module") -> temp_pblue1


png(filename = glue("results_figures/MEblue_CORT_{my_tissue}_biggerfont.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
temp_pblue1
invisible(dev.off())

ME_df %>% 
  ggplot(aes(status,MEblue, color = status, fill = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 jitter.color = NA,
                 alpha = 0.4,
                 jitter.size = 2.2,
                 jitter.height = 0.02, jitter.width = 0.1, 
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",
        plot.title = element_text(size = 7)) +
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Liver-Blue module") -> temp_pblue2


png(filename = glue("results_figures/MEblue_status_{my_tissue}.png"),
    width = standard_size, height = standard_size, units = "cm", res = 600)
temp_pblue2
invisible(dev.off())



ME_df %>% 
  ggplot(aes(cort_post,MEblue,color = status, fill = status))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(shape = 21, alpha = 0.6, size = 2.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(size = 10),
        # legend.position = c(0.2,0.83),
        legend.position = "none",)+
  labs(x = "Plasma CORT (GD14)",
       y = "Module eigengene",
       title = "Liver-Blue module",
       fill ="",
       color = "") -> temp_pblue3


png(filename = glue("results_figures/MEblue_INT2_{my_tissue}_bf.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
temp_pblue3
invisible(dev.off())

ME_df %>% 
  ggplot(aes(MLR_GD14,MEblue,color = status, fill = status))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(shape = 21, alpha = 0.6, size = 2.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 7),
        # legend.position = c(0.2,0.83),
        legend.position = "none",
        legend.key.size = unit(0.55, 'cm'),
        legend.text = element_text(size = 5))+
  labs(x = "Monocyte-lymphocyte ratio (GD14)",
       y = "Module eigengene",
       title = "Liver-Blue module",
       fill ="",
       color = "") -> temp_pblue4


png(filename = glue("results_figures/MEblue_MLR_{my_tissue}.png"),
    width = standard_size, height = standard_size, units = "cm", res = 600)
temp_pblue4
invisible(dev.off())

# ME_df %>%
  # ggplot(aes(MLR_GD14,cort_post,color = status))+geom_point()


# midnight blue
ME_df %>% 
  ggplot(aes(Mono_GD14,MEmidnightblue))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(size = 1.6, shape = 21, alpha = 0.6,size = 2.2)+
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 7))+
  labs(x = "Monocyte (%) (GD14)",
       y = "Module eigengene",
       title = "Liver-Midnightlue module") -> temp_pblue1


png(filename = glue("results_figures/MEmidnightblue_mono_{my_tissue}.png"),
    width = standard_size, height = standard_size, units = "cm", res = 600)
temp_pblue1
invisible(dev.off())


# black
ME_df %>% 
  ggplot(aes(status,MEblack, color = status, fill = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 jitter.color = NA,
                 alpha = 0.4,
                 jitter.size = 2.2,
                 jitter.height = 0.02, jitter.width = 0.1, 
                 errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",
        plot.title = element_text(size = 7)) +
  labs(x = "Social status",
       y = "Module eigengene",
       title = "Liver-Black module") -> temp_pblue2


png(filename = glue("results_figures/MEblack_status_{my_tissue}.png"),
    width = standard_size, height = standard_size, units = "cm", res = 600)
temp_pblue2
invisible(dev.off())



ME_df %>% 
  ggplot(aes(cort_post,MEblack,color = status))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8, se = F)+
  geom_point(size = 1.6, shape = 21, alpha = 0.6, size = 2.2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 7),
        legend.position = c(0.2,0.83),
        legend.key.size = unit(0.55, 'cm'),
        legend.text = element_text(size = 5))+
  labs(x = "Plasma CORT (GD14)",
       y = "Module eigengene",
       title = "Liver-Black module",
       fill ="",
       color = "") -> temp_pblue3


png(filename = glue("results_figures/MEblack_INT_{my_tissue}.png"),
    width = standard_size, height = standard_size, units = "cm", res = 600)
temp_pblue3
invisible(dev.off())


# Spleen - brown and green ME ~ CORT =================================

# make sure to rerun the top codes for spleen dataset before you run these below 
ME_df %>% 
  ggplot(aes(cort_post,MEbrown))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(size = 1.6, shape = 21, alpha = 0.6, size = 2.2)+
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(size = 10))+
  labs(x = "Plasma CORT (GD14)",
       y = "Module eigengene",
       title = "Spleen-Orange module") -> temp_p4


png(filename = glue("results_figures/MEbrown_{my_tissue}_biggerfont.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
temp_p4
invisible(dev.off())


ME_df %>% 
  ggplot(aes(cort_post,MEgreen))+
  geom_smooth(method = "lm", alpha = 0.2, size = 0.8)+
  geom_point(size = 1.6, shape = 21, alpha = 0.6, size = 2.2)+
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(size = 10))+
  labs(x = "Plasma CORT (GD14)",
       y = "Module eigengene",
       title = "Spleen-Darkgreen module") -> temp_p5


png(filename = glue("results_figures/MEgreen_{my_tissue}_biggerfont.png"),
    width = standard_size+0.7, height = standard_size, units = "cm", res = 600)
print(temp_p5)
invisible(dev.off())

























moduleColors %>% 
  table() %>% 
  as.data.frame() %>% arrange(Freq)  -> modnum 
colnames(modnum) <- c("module","count")
str(modnum)

modnum$module <- factor(modnum$module,levels=rev(unique(modnum$module[order(modnum$count)])))
mycol <-modnum$module %>% as.character() # for some reason this has to be character again to be ordered properly in the figure...!! 


modnum %>%
  ggplot(aes(y = module, x = count, fill = module))+
  geom_bar(stat = 'identity', color = NA)+
  geom_text(aes(label = count),hjust = 1.1, size = 1.4, angle = 180)+
  xlim(0, max(modnum$count)+75)+
  scale_fill_manual(values = rev(mycol))+
  theme_minimal(base_size = 7)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 180, margin = margin(r=-6)),
        panel.grid = element_blank())+
  labs(x = "",
       y = "") -> temp_p
temp_p

png(filename = glue("results_figures/modulesize_{my_tissue}_upsidedown.png"),
    width = 6, height = my_height, units = "cm", res = 600)
temp_p
invisible(dev.off())


# without grey module 



my_tissue = "Liver"
my_tissue = "Spleen"

limma_list<- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS")) 


png(filename = glue("results_figures/EnVol_{my_tissue}_status.png"),
        width = 12, height = 11, units = "cm", res = 600)
EnhancedVolcano(limma_list$status,
                lab = limma_list$status$symbol,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue}: DEG by social status"),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                shape = 21,
                pointSize = 2.5,
                labSize = 2.0)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.5),breaks = c(0,1,2,3,4))+
  theme_bw(base_size = 7)+
  labs(color = "")+
  theme(legend.position = "top")
invisible(dev.off())

png(filename = glue("results_figures/EnVol_{my_tissue}_cort.png"),
    width = 12, height = 11, units = "cm", res = 600)
EnhancedVolcano(limma_list$cort,
                lab = limma_list$cort$symbol,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue}: DEG by CORT"),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                shape = 21,
                pointSize = 2.5,
                labSize = 2.0)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.5),breaks = c(0,1,2,3,4))+
  theme_bw(base_size = 7)+
  labs(color = "")+
  theme(legend.position = "top")
invisible(dev.off())

# pval and no label 

keyvals <- ifelse(
  limma_list$status$logFC < 0 & limma_list$status$P.Value<0.05, 'orange',
  ifelse(limma_list$status$logFC > 0 & limma_list$status$P.Value<0.05, 'purple4',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'

png(filename = glue("results_figures/EnVol_{my_tissue}_status2.png"),
    width = 12, height = 11, units = "cm", res = 600)

EnhancedVolcano(limma_list$status,
                lab = limma_list$status$symbol,
                selectLab=NA,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue}: DEG by social status"),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "blank",
                labSize = 0.0)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.5),breaks = c(0,1,2,3,4))+
  theme_bw(base_size = 7)+
  labs(color = "")+
  theme(legend.position = "none")

invisible(dev.off())



keyvals <- ifelse(
  limma_list$cort$logFC < 0 & limma_list$cort$P.Value<0.05, 'orange',
  ifelse(limma_list$cort$logFC > 0 & limma_list$cort$P.Value<0.05, 'purple4',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'purple4'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'orange'] <- 'low'

png(filename = glue("results_figures/EnVol_{my_tissue}_cort2.png"),
    width = 12, height = 11, units = "cm", res = 600)

EnhancedVolcano(limma_list$cort,
                lab = limma_list$cort$symbol,
                selectLab=NA,
                x = 'logFC',
                y = 'P.Value',
                title = glue("{my_tissue}: DEG by CORT"),
                pCutoff = 0.05,
                FCcutoff = 0.0,
                cutoffLineType = 'blank',
                hline = c(0.05),
                hlineCol = c('grey'),
                hlineType = c( 'dotted'),
                hlineWidth = c(1),
                colCustom = keyvals,
                shape = 21,
                pointSize = 2.5,
                labCol = "blank",
                labSize = 0.0)+
  scale_x_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  scale_y_continuous(limits = c(-0.1,4.5),breaks = c(0,1,2,3,4))+
  theme_bw(base_size = 7)+
  labs(color = "")+
  theme(legend.position = "none")

invisible(dev.off())

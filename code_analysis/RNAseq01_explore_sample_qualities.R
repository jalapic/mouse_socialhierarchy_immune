colSums(rnaseq_rawcounts[,-1]) %>% 
  as_tibble_row() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'sample_id') %>% 
  rename(genecounts = V1) -> genecounts

genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  theme_base() 
# ggsave(p_genecounts,filename = 'results_figures/p_genecounts.png',width = 8, height = 5)

genecounts %>% 
  filter(genecounts < 800000) 

genecounts %>% 
  filter(genecounts > 3000000)

genecounts %>% 
  ggplot(aes(sample_id,genecounts))+
  geom_bar(stat = 'identity')+
  coord_flip()

genecounts %>%
  mutate(sampleID = gsub(".trim.sam.counts","",sample_id)) %>% 
  left_join(rnaseq_sampleid) %>% 
  left_join(alldata %>% select(subjectID, domgroup)) %>% 
  filter(domgroup!="Subdominant") %>% 
  ggplot(aes(domgroup,genecounts))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 0.3,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.6)) +
  facet_wrap(~tissue, scales = "free_y")



my_tissue = 'Liver'
my_tissue = 'Spleen'
# filter_counts = 100
# filter_counts = 50 # for WGCNA dds 
filter_counts = 20 # for PCA
# filter_counts = 10

coldata = rnaseq_sampleid %>% 
  filter(subjectID!='G.5') %>% 
  filter(tissue == my_tissue) %>% 
  mutate(subjectID = as.character(subjectID)) %>% 
  left_join (.,all_behavior %>% 
               dplyr::select(subjectID, domgroup), by = 'subjectID') %>% 
  column_to_rownames(var = 'sampleID') %>% 
  mutate(sampleID = row.names(.))


colnames(rnaseq_rawcounts) <- gsub(".trim.sam.counts","",colnames(rnaseq_rawcounts))

# tidy count data format  
df <- rnaseq_rawcounts[,colnames(rnaseq_rawcounts) %in% coldata$sampleID]
countData <- as.matrix(df) 
dim(countData)
rownames(countData) <- rnaseq_rawcounts$ensgene
countData <- countData[rowSums(countData > filter_counts) > round((length(df))*0.9), ]

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~domgroup)

dds <- DESeq(dds)
dim(dds)
# saveRDS(dds, glue("results_RNAseqRDS/dds_{my_tissue}.RDS"))


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = c("domgroup"), returnData = T) %>%
  left_join(coldata %>% 
              rename(name = sampleID))-> d

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("purple4", "orange")) + 
  labs(x = paste0("PC1: ", 
              round(attributes(d)$percentVar[1] * 100), 
              "% variance"),
       y = paste0("PC2: ", 
              round(attributes(d)$percentVar[2] * 100), 
              "% variance"),
       color = "Social status",
       fill = "Social status") +
  theme(legend.position = c(0.76,0.85),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11))+
  ggtitle(glue('{my_tissue} PCA plot')) -> p



theme_set(theme_bw(base_size = 12.5))

png(filename = glue("results_figures/PCA_{my_tissue}_biggerfont.png"),
    width = 8.5, height = 9, units = "cm", res = 600)
p
invisible(dev.off())



ggplot(data = d %>% 
         left_join(alldata %>% select(subjectID, cort_post)),
       aes_string(x = "PC1", y = "PC2", color = "cort_post")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_viridis(option = "D") +
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "CORT (ng/ml)",
       fill = "CORT (ng/ml)") +
  theme(legend.position = c(0.88,0.83),
        legend.key.size = unit(0.35, 'cm'),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))+
  ggtitle(glue('{my_tissue} PCA plot')) -> p2

p2

theme_set(theme_bw(base_size = 10))

png(filename = glue("results_figures/PCA_cort_{my_tissue}.png"),
    width = 8.5, height = 9, units = "cm", res = 600)
p2
invisible(dev.off())





# pca for both status and cort  ================

ggplot(data = d %>% 
         left_join(alldata %>% select(subjectID, cort_post)),
       aes_string(x = "PC1", y = "PC2", shape = 'group',color = "cort_post", fill = "cort_post")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_viridis(option = "D") +
  scale_fill_viridis(option = "D") +
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       shape ="Social status",
       color = "CORT (ng/ml)",
       fill = "CORT (ng/ml)") +
  theme(legend.key.size = unit(0.25, 'cm'),
      legend.margin = margin(-0.1,0,0,0, unit="cm"),
      # legend.position = c(0.88,0.73),
        legend.text = element_text(size=4),
        legend.title = element_text(size=5))+
  # guides(shape = guide_legend(override.aes = list(size = 2)))+
  ggtitle(glue('{my_tissue} PCA plot')) -> p3

p3

theme_set(theme_bw(base_size = 10))

png(filename = glue("results_figures/PCA_both_{my_tissue}.png"),
    width = 10.8, height = 9, units = "cm", res = 600)
p3
invisible(dev.off())

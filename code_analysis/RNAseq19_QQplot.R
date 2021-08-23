# Got this idea from Jenny Tung's S4SN 2020* talk! 

my_tissue = "Liver"
my_tissue = "Spleen"

if(my_tissue == "Liver"){limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))
} else {limma_list <- readRDS(glue("results_RNAseqRDS/limma_{my_tissue}.RDS"))}

limma_list %>% 
  map(~arrange(., P.Value)) %>% 
  map(~head(.,20))

p <- limma_list$status$P.Value

nn = length(p)
xx =  -log10((1:nn)/(nn+1))

p_cort <- limma_list$cort$P.Value

length(p) == length(p_cort) # TRUE
sum(p == p_cort)
df <- cbind(exp = xx, 
            Status = -sort(log10(p)),
            CORT = -sort(log10(p_cort))) %>% 
  as.data.frame() %>% 
  gather(key,value,2:3) %>% 
  mutate(key = factor(key, levels = c("Status","CORT")))

lim = log10(5000)


png(filename = glue("results_figures/QQplot_{my_tissue}_biggerfont.png"),
    width = 8, height = 9, units = "cm", res = 600)
df %>%
  filter(value < lim) %>% 
  ggplot(aes(exp, value, color = key))+
  geom_point(shape  = 21,size =1)+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype ="dashed")+
  theme_bw(base_size = 12.5)+
  labs(x = expression(-log[10]~p-values~(expected)),
       y = expression(-log[10]~p-values~(observed)),
       title = glue("{my_tissue}"),
       color = "",
       fill = "")+
  theme(legend.position = c(0.8,0.2),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3.3)))+
  scale_color_manual(values = c("#E7B800","#00AFBB" ))+
  scale_fill_manual(values = c("#E7B800","#00AFBB" ))+
  scale_x_continuous(limits = c(0, 3.5), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(0, 3.5), expand=expansion(mult=c(0,0.0)))

dev.off()
  

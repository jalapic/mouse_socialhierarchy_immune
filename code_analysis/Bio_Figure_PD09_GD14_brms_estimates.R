# 
# saveRDS(plotdf_pd09_pairstatus,"results_statRDS/plotdf_pd09_pairstatus.RDS")
# saveRDS(plotdf_pd09_predict,"results_statRDS/plotdf_pd09_predict.RDS")
# saveRDS(plotdf_pd09_predict_glick,"results_statRDS/plotdf_pd09_predict_glick.RDS")
# saveRDS(plotdf_gd14_ds,"results_statRDS/plotdf_gd14_ds.RDS")
# saveRDS(plotdf_gd14_moglick,"results_statRDS/plotdf_gd14_moglick.RDS")
# saveRDS(plotdf_gd14_glickorating,"results_statRDS/plotdf_gd14_glickorating.RDS")
# saveRDS(plotdf_gd14_domcert,"results_statRDS/plotdf_gd14_domcert.RDS")
options(backup_options)

plotdf_pd09_pairstatus <- readRDS("results_statRDS/plotdf_pd09_pairstatus.RDS")
plotdf_pd09_predict <- readRDS("results_statRDS/plotdf_pd09_predict.RDS")
# plotdf_pd09_predict_glick <- readRDS("results_statRDS/plotdf_pd09_predict_glick.RDS")
plotdf_gd14_ds <- readRDS("results_statRDS/plotdf_gd14_ds.RDS")
plotdf_gd14_moglick <- readRDS("results_statRDS/plotdf_gd14_moglick.RDS")
plotdf_gd14_glickorating <- readRDS("results_statRDS/plotdf_gd14_glickorating.RDS")
# plotdf_gd14_domcert <- readRDS("results_statRDS/plotdf_gd14_domcert.RDS")

plotdf_pd09_pairstatus %>% 
  mutate(b_pair_statussub =(-1)*b_pair_statussub,
         .upperx = .upper,
         .upper = (-1)*.lower,
         .lower = (-1)*.upperx) -> plotdf_pd09_pairstatus


mycellname <- data.frame (
  cell_type = plotdf_pd09_pairstatus %>% select(cell_type) %>% unique(),
  pretty_cell_type = c("B cells", "CD4/CD8", "Cytotoxic T", "DCs",
                       "Helper T", "Lymphoid DCs", "Macrophages","MLR",
                       "Monocytes", "Myeloid DCs","Neutrophils","NK cells",
                       "NLR","T cells")) %>% 
  mutate(pretty_cell_type = factor(pretty_cell_type,
                                   levels = rev(c("Macrophages","Monocytes","Neutrophils",
                                                  "NK cells",
                                                  "Lymphoid DCs","Myeloid DCs", "DCs",
                                                  "B cells","T cells","Helper T","Cytotoxic T",
                                                  "CD4/CD8","MLR","NLR"))))






# plot 1 
plotdf_pd09_pairstatus[,1:3] <- (-1)*plotdf_pd09_pairstatus[,1:3]
colnames(plotdf_pd09_pairstatus)[2:3] <- c( ".upper", ".lower")
plotdf_pd09_pairstatus %>% 
  mutate(b_pair_statussub = (-1)* b_pair_statussub) %>%
  mutate(`.upperx` = (-1)*`.lower`,
         `.lowerx` = (-1)*`.upper`) %>%
  mutate(`.upper` = `.upperx`,
         `.lower` = `.lowerx`) -> plotdf_pd09_pairstatus
plotdf_pd09_pairstatus %>%
mutate(myx = b_pair_statussub) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_pd09_pairstatus) %>% 
  mutate(myx = b_pair_statussub) %>% #
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs")%>% 
  filter(pretty_cell_type != "CD4/CD8")  -> plotdf1

plotdf1 %>% 
  ggplot(aes(y = pretty_cell_type, x = myx ,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf1 %>% filter(.width == 0.66),
            aes(y = pretty_cell_type, yend = pretty_cell_type,
                x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf1 %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Subordinate       --       Estimates (median ± 95% CI)       --       Dominant",
      y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf1 %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -1.65), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf1 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -1.65), hjust = "inward", color = "black", )+ #
  xlim(c(-1.65,0.78))+
  ggtitle("Effect of dominant-subordinate status in pairs") -> plot1

print(plot1)

# 
# plotdf1 %>% 
#   ggplot(aes(y = pretty_cell_type, x = myx,color = sig, fill = sig))+ 
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
#   geom_pointintervalh()+
#   theme_base() +
#   labs(x = "Estimates (median ± 95% CI)",
#        y = "")+
#   theme(legend.position = "none",
#         axis.title.x = element_text(size = 12),
#         plot.title = element_text(size = 14))+
#   scale_color_manual(values = c("black","red"))+
#   geom_text(data = plotdf1 %>% filter(sig == "Sig"), 
#             aes(label = statsummary, x= -1.35), hjust = "inward", color = "red", fontface = "bold")+ #
#   geom_text(data = plotdf1 %>% filter(sig != "Sig"), 
#             aes(label = statsummary, x= -1.35), hjust = "inward", color = "black", )+ #
#   xlim(c(-1.35,1.05))+
#   ggtitle("Effect of dominant-subordinate status in pairs (Dom - Sub)")


# plot 2 
plotdf_pd09_predict %>% 
  mutate(myx = b_perc) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_pd09_predict) %>% 
  mutate(myx = b_perc) %>% #
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs")%>% 
  filter(pretty_cell_type != "CD4/CD8")  -> plotdf2

plotdf2 %>% 
  ggplot(aes(y = pretty_cell_type, x = myx,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf2 %>% filter(.width == 0.66),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf2 %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Estimates (median ± 95% CI)",
       y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf2 %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -1.8), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf2 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -1.8), hjust = "inward", color = "black", )+ #
  xlim(c(-1.8,0.9))+
  ggtitle("Predictabiliy of final rank by immunophenotypes on PD09") -> plot2

print(plot2)



# plot 3
plotdf_gd14_ds %>% 
  mutate(myx = b_ds) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_gd14_ds) %>% #
  mutate(myx = b_ds) %>% 
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs") %>% 
  filter(pretty_cell_type != "CD4/CD8")-> plotdf3

plotdf3 %>% 
  ggplot(aes(y = pretty_cell_type, x = myx,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf3 %>% filter(.width == 0.66),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf3 %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Subordinate       --       Estimates (median ± 95% CI)       --       Dominant",
       y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf3 %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -0.32), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf3 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -0.32), hjust = "inward", color = "black", )+ #
  xlim(c(-0.32,0.23))+
  ggtitle("Immunophenotypes on GD14 ~ final David's score") -> plot3


print(plot3)



# Let's print
theme_set(theme_bw(base_size = 10))

png(filename = "results_figures/flow_pairstatus.png",
    width = 18, height = 12, units = "cm", res = 600)
print(plot1)
invisible(dev.off())

png(filename = "results_figures/flow_predict_finalrank.png",
    width = 18, height = 12, units = "cm", res = 600)
print(plot2)
invisible(dev.off())

png(filename = "results_figures/flow_group_davidscore.png",
    width = 18, height = 12, units = "cm", res = 600)
print(plot3)
invisible(dev.off())









# plot 4

plotdf_gd14_domcert %>% 
  mutate(myx = b_ind_cert) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_gd14_domcert) %>% #
  mutate(myx = b_ind_cert) %>% 
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs")  -> plotdf4


plotdf4 %>% 
  ggplot(aes(y = pretty_cell_type, x = myx,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf4 %>% filter(.width == 0.66),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf4 %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Uncertain       --       Estimates (median ± 95% CI)       --       Certain",
       y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf4 %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -14), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf4 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -14), hjust = "inward", color = "black", )+ #
  xlim(c(-14,12))+
  ggtitle("Immunophenotypes on GD14 ~ Dominance certainty") -> plot4

print(plot4)



# mywidth = 7.2
# myheight = 5.8
# ggsave(plot1, filename = "results_figures/flow_brm_pairdomsub.png",width = mywidth,height = myheight)  
# ggsave(plot2, filename = "results_figures/flow_brm_predictfinalrank.png",width = mywidth,height = myheight)  
# ggsave(plot3, filename = "results_figures/flow_brm_davidscore.png",width = mywidth,height = myheight)  
# ggsave(plot4, filename = "results_figures/flow_brm_plot4.png",width = mywidth,height = myheight)  
# 
# x <- plot_grid(plot1, plot2, plot3, plot4, labels = c('A', 'B','C','D'), label_size = 15, nrow =1)
# xx <- plot_grid(plot1, plot2, plot3, plot4, labels = c('A', 'B','C','D'), label_size = 15, nrow =2)
# 
# ggsave(x, filename = "results_figures/flow_brm_all.png",width = mywidth*4,height = myheight)  
# ggsave(xx, filename = "results_figures/flow_brm_all2.png",width = mywidth*2,height = myheight*2)  










# GD14 flow data - for rank, multiply ?  
plotdf_gd14_moglick %>% 
  mutate(myx = bsp_moglicko_rank) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_gd14_moglick) %>% #
  mutate(myx = bsp_moglicko_rank) %>% 
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs")  -> plotdf_moglick

plotdf_moglick %>% 
  ggplot(aes(y = pretty_cell_type, x = myx,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf_moglick %>% filter(.width == 0.66),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf_moglick %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Subordinate       --       Estimates (median ± 95% CI)       --       Dominant",
       y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf_moglick %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -0.32), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf_moglick %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -0.32), hjust = "inward", color = "black", )+ #
  xlim(c(-0.32,0.23))+
  ggtitle("Immunophenotypes on GD14 ~ final social rank") -> plot_moglick


print(plot_moglick)
ggsave(plot_moglick, filename = "results_figures/flow_brm_GD14_moglick.png",width = mywidth,height = myheight)  






# GD14 ~ Glicko ratings
plotdf_gd14_glickorating %>% 
  mutate(myx = b_scaled_Rating) %>% 
  filter(`.width` == 0.95) %>% 
  mutate(statsummary = glue ('{round(myx,2)} [{round(.lower,2)}, {round(.upper,2)}]')) %>% 
  mutate(sig = ifelse(.lower*.upper >= 0, "Sig", "")) %>% 
  select(cell_type,sig,statsummary) %>% 
  right_join(plotdf_gd14_glickorating) %>% #
  mutate(myx = b_scaled_Rating) %>% 
  left_join(mycellname) %>% 
  filter(pretty_cell_type != "DCs")  -> plotdf5

plotdf5 %>% 
  ggplot(aes(y = pretty_cell_type, x = myx,color = sig))+ 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")+
  geom_point(size = 3.5)+
  geom_segment(data = plotdf5 %>% filter(.width == 0.66),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 1.5)+
  geom_segment(data = plotdf5 %>% filter(.width == 0.95),
               aes(y = pretty_cell_type, yend = pretty_cell_type,
                   x = .lower, xend = .upper), size = 0.75)+
  theme_base() +
  labs(x = "Subordinate       --       Estimates (median ± 95% CI)       --       Dominant",
       y = "")+
  theme(legend.position = "none",
        plot.background = element_blank(),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14))+
  scale_color_manual(values = c("black","red"))+
  geom_text(data = plotdf5 %>% filter(sig == "Sig"), 
            aes(label = statsummary, x= -0.62), hjust = "inward", color = "red", fontface = "bold")+ #
  geom_text(data = plotdf5 %>% filter(sig != "Sig"), 
            aes(label = statsummary, x= -0.62), hjust = "inward", color = "black", )+ #
  # xlim(c(-0.32,0.23))+
  ggtitle("Immunophenotypes on GD14 ~ Glicko rating (scaled)") -> plot5


print(plot5)

mywidth = 7.8
myheight = 8
ggsave(plot5, filename = "results_figures/flow_brm_scaled_glicko.png",width = mywidth,height = myheight)  


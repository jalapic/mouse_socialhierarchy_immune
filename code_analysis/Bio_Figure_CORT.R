
source('code_carpentry/00_preset_ggplot_options.R')

cortdf <- readRDS("results_statRDS/cortdf.RDS")


# Pair housing Dom vs. Sub

cortdf %>% 
  mutate(buddyID = paste(batch, mouseID, sep = "-")) %>% 
  ggplot(aes(pair_status, cort_pre, group = buddyID, color = pair_status, fill = pair_status))+
  geom_line(alpha = 0.2, size = 0.8, color = 'grey15')+
  geom_point(shape = 21, alpha = 0.3,size = 2)+
  scale_x_discrete(labels=c("Dom","Sub"))+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme(legend.position = "none") +
  labs(x = "Pair housing social status", y = "Corticosterone level (ng/ml)") -> plot_paircort

print(plot_paircort)



# Group housing 
cortdf %>% 
  ggplot(aes(domgroup,cort_post, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # scale_color_brewer(palette = "Dark2")+
  # scale_fill_brewer(palette = "Dark2")+
  scale_color_manual(values = c("purple4","#21908CFF","orange"))+
  scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
  theme_bw(base_size = 11.2)+
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.02))) +
  labs(x = "Social status", y = "Corticosterone level (ng/ml)")+
  # ylim(0,680)+
  facet_wrap(~despotism_binary) -> plot_groupcort

print(plot_groupcort) # awesome, pretty much replicating the hormone paper



png(filename = "results_figures/grouphousing_cort_biggerstripfont.png",
    width = 14, height = 10, units = "cm", res = 600)
print(plot_groupcort)
invisible(dev.off())




# Pre - Post change by social status 
cortdf %>% 
  gather(key,value, c('cort_pre','cort_post')) %>% 
  mutate(key = ifelse(key == "cort_post", "Post","Pre")) %>% 
  mutate(key = factor(key, levels = c('Pre','Post'))) %>% 
  mutate(pair_status = ifelse(pair_status == "dom", "Dominant","Subordinate")) %>% 
  ggplot(aes(key,value,group = subjectID, color = pair_status, fill = pair_status))+
  geom_line(alpha = 0.2, size = 1)+
  geom_point(shape = 21, alpha = 0.3,size = 2)+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  facet_wrap(~domgroup)+
  theme(legend.position = "right")+
  labs(x = "", y = "Corticosterone level (ng/ml)",
       color = 'Pair housing\nsocial status',
       fill = 'Pair housing\nsocial status') -> plot_prepostcort

print(plot_prepostcort)



# Facet by pair housing social

cortdf %>% 
  ggplot(aes(domgroup,cort_post, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none") +
  labs(x = "Social status", y = "Corticosterone level (ng/ml)")+
  facet_wrap(~pair_status, 
             labeller = labeller(pair_status = 
                                   c(dom = "Pair housing status: Dominant",
                                     sub = "Pair housing status: Subordinate"))) -> plot_groupcort_pairstatus

print(plot_groupcort_pairstatus)

# Print

theme_set(theme_bw(base_size = 10))

png(filename = "results_figures/pairhousing_cort.png",
    width = 7, height = 10, units = "cm", res = 600)
print(plot_paircort)
invisible(dev.off())


png(filename = "results_figures/grouphousing_cort.png",
    width = 14, height = 10, units = "cm", res = 600)
print(plot_groupcort)
invisible(dev.off())


png(filename = "results_figures/prepost_cort.png",
    width = 14, height = 10, units = "cm", res = 600)
print(plot_prepostcort)
invisible(dev.off())


png(filename = "results_figures/grouphousing_cort_pairstatus.png",
    width = 14, height = 10, units = "cm", res = 600)
print(plot_groupcort_pairstatus)
invisible(dev.off())




# # Pair housing Dom vs. Sub
cort_pre <- brm(cort_pre ~ pair_status + (1|pairID),
                data = cortdf %>% 
                  mutate(pairID = paste(batch,mouseID,sep = "_")),
                prior = set_prior('normal(0, 30)'),
                file = "results_statRDS/cort_pre.Rds") # should be already saved (Bio11_cort.R)

conditional_effects(cort_pre)$pair_status %>% 
  rename(estimate = estimate__,
         lower = lower__,
         upper = upper__) -> brms_temp
cortdf %>% 
  mutate(buddyID = paste(batch, mouseID, sep = "-")) %>% 
  ggplot(aes(pair_status, cort_pre, group = buddyID, color = pair_status, fill = pair_status))+
  geom_line(alpha = 0.2, size = 0.8, color = 'grey15')+
  geom_point(shape = 21, alpha = 0.3,size = 2)+
  annotate("pointrange",
           x = 2.2 ,
           y = brms_temp %>% filter(pair_status != "dom") %>% .$estimate, 
           ymin = brms_temp %>% filter(pair_status != "dom") %>% .$lower, 
           ymax = brms_temp %>% filter(pair_status != "dom") %>% .$upper,
           size = 0.7,
           alpha=0.9)+
  annotate("pointrange",
           x = 0.8, 
           y = brms_temp %>% filter(pair_status == "dom") %>% .$estimate, 
           ymin = brms_temp %>% filter(pair_status == "dom") %>% .$lower, 
           ymax = brms_temp %>% filter(pair_status == "dom") %>% .$upper,
           size = 0.7,
           alpha=0.9) +
  scale_x_discrete(labels=c("Dom","Sub"))+
  scale_color_manual(values = c("purple4","orange"))+
  scale_fill_manual(values = c("purple4","orange"))+
  theme(legend.position = "none") +
  ylim(0,395)+
  labs(x = "Pair housing social status", y = "Corticosterone level (ng/ml)") -> plot_paircort

print(plot_paircort)


theme_set(theme_bw(base_size = 10))

png(filename = "results_figures/pairhousing_cort_new.png",
    width = 7, height = 10, units = "cm", res = 600)
print(plot_paircort)
invisible(dev.off())







# Group housing 
cortdf %>% 
  ggplot(aes(domgroup,cort_post, color = domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.size = 2,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # scale_color_brewer(palette = "Dark2")+
  # scale_fill_brewer(palette = "Dark2")+
  scale_color_manual(values = c("purple4","#21908CFF","orange"))+
  scale_fill_manual(values = c("purple4","#21908CFF","orange"))+
  theme(legend.position = "none") +
  labs(x = "", y = "Corticosterone level (ng/ml)") -> plot_groupcortx

print(plot_groupcortx) # awesome, pretty much replicating the hormone paper



png(filename = "results_figures/group_cort_altogether.png",
    width = 8, height = 8, units = "cm", res = 600)
print(plot_groupcortx)
invisible(dev.off())

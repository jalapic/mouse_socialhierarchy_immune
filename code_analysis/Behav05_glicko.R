
# Get glicko scores
df.glickos <-  lapply(df.groupxx, get_glickos, cval=3)


# Make glicko plot
pglickos<-NULL
for(i in 1:length(df.glickos)){
  pglickos[[i]]<- plotglicko(df.glickos[[i]], ylim1=1300,ylim2=3300,  linewd=0.8)+
    geom_hline(yintercept = 2200, linetype = "dashed", color = "grey50" )+
    theme_classic()+
    theme(legend.position = "none")+
    ggtitle(glue("Cohort {LETTERS[[i]]}"))+
    theme(axis.text = element_text(color="#3C3C3C", size=rel(1.1)),
          axis.title = element_text(color="#3C3C3C", size=rel(0.9)),
          title = element_text(size=rel(1.3)),
          legend.position = "none"   )
  
}

# plotglickos <- grid.arrange(pglickos[[1]],pglickos[[2]],pglickos[[3]],pglickos[[4]],
#                             pglickos[[5]],pglickos[[6]],pglickos[[7]],pglickos[[8]],
#                             pglickos[[9]],pglickos[[10]],pglickos[[11]],pglickos[[12]],
#                             nrow=3 )
png(filename = "results_figures/Glicko_all_biggerfont.png",
    width = 27, height = 21, units = "cm", res = 600)
grid.arrange(pglickos[[1]],pglickos[[2]],pglickos[[3]],pglickos[[4]],
             pglickos[[5]],pglickos[[6]],pglickos[[7]],pglickos[[8]],
             pglickos[[9]],pglickos[[10]],pglickos[[11]],pglickos[[12]],
             nrow=3 )

invisible(dev.off())


# individual glicko graph ================================================

i = 9

png(filename = "results_figures/Glicko_cohortI_biggerfont.png",
    width = 12, height = 12, units = "cm", res = 600)
plotglicko(df.glickos[[i]],
           ylim1=1300,ylim2=3300,  linewd=0.7)+
  geom_hline(yintercept = 2200, linetype = "dashed", color = "grey50" )+
  theme_classic()+
  labs(x = "Aggressive interaction events",
      y = "Dominance score (Glicko rating)")+
  theme(axis.text = element_text(color="#3C3C3C", size=rel(1.35)),
        axis.title = element_text(color="#3C3C3C", size=rel(1.25)),
        title = element_text(size=rel(1.3)),
        legend.position = "none"   )
invisible(dev.off())



i = 7 # displacement example 
png(filename = "results_figures/Glicko_cohortG.png",
    width = 12, height = 12, units = "cm", res = 600)
plotglicko(df.glickos[[i]],
           ylim1=1300,ylim2=3300,  linewd=0.7)+
  geom_hline(yintercept = 2200, linetype = "dashed", color = "grey50" )+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = "Aggressive interaction events",
       y = "Dominance score (Glicko rating)")
invisible(dev.off())

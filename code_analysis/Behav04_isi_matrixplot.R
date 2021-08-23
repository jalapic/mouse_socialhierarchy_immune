### Matrix Plot.

# Organize matrices by David's Scores
m.wlds <- lapply(m.wl, compete::org_matrix, method="ds")


# Organize matrices by ISI
m.wlds <- lapply(m.wl, compete::org_matrix, method="ds")
m.dids <- lapply(m.di, compete::org_matrix, method="ds")
m.isi <- lapply(m.wlds, function(x) compete::isi98(x)$best_order)


bimat<-socio<-vector('list',12)




theme_get()


for(i in 1:12){
  mm <- m.wlds[[i]][m.isi[[i]],m.isi[[i]]]
  rownames(mm)<-colnames(mm)<-1:10
  socio[[i]] <- matrixplot(mm)+ ggtitle(glue("Cohort {LETTERS[i]}"))  +
    ylab("Winner's Rank") +
    xlab("Loser's Rank")
 }

# 
# for(i in 1:12){
#   mm <- m.wlds[[i]][m.isi[[i]],m.isi[[i]]]
#   rownames(mm)<-colnames(mm)<-1:10
#   bimat[[i]] <- matrixplot0(mm) + ggtitle(LETTERS[i])  +
#     ylab("Winner's Rank") +
#     xlab("Loser's Rank")
# }
# 
# matrices=grid.arrange(socio[[1]],socio[[2]],socio[[3]],socio[[4]],socio[[5]],socio[[6]],
#                       socio[[7]],socio[[8]],socio[[9]],socio[[10]],socio[[11]],socio[[12]],
#                       bimat[[1]],bimat[[2]],bimat[[3]],bimat[[4]],bimat[[5]],bimat[[6]],
#                       bimat[[7]],bimat[[8]],bimat[[9]],bimat[[10]],bimat[[11]],bimat[[12]],
#                       nrow=2)

sociomatrices=grid.arrange(socio[[1]],socio[[2]],socio[[3]],socio[[4]],socio[[5]],socio[[6]],
                      socio[[7]],socio[[8]],socio[[9]],socio[[10]],socio[[11]],socio[[12]],
                      ncol=4)


ggsave("results_figures/sociomatrices_biggerfont.png", sociomatrices, width=15,height=12)



# bimatrices=grid.arrange(bimat[[1]],bimat[[2]],bimat[[3]],bimat[[4]],bimat[[5]],bimat[[6]],
#                       bimat[[7]],bimat[[8]],bimat[[9]],bimat[[10]],bimat[[11]],bimat[[12]],
#                       ncol=4)
# 
# ggsave("results_figures/bimatrices.png", bimatrices, width=15,height=12)


# individual sociomatrix ==========================================
i = 9

png(filename = "results_figures/sociomat_cohortI_new_biggerfont.png",
    width = 12, height = 12, units = "cm", res = 600)

socio[[i]]+
    labs(x = "Loser's rank",
         y = "Winner's rank")+
    theme(plot.title = element_blank(), 
          panel.background = element_blank(), 
          plot.background = element_blank(), 
          panel.grid.major.y = element_line(color = "gray75",linetype = 'dotted'), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          text = element_text(color="gray20", size=10),
          axis.text = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="gray20", size=rel(1.45)),
          axis.text.y = element_text(color="gray20", size=rel(1.45)),
          axis.title.x = element_text(size=rel(1.2), vjust=0),
          axis.title.y = element_text(size=rel(1.2), vjust=1),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
          
invisible(dev.off())

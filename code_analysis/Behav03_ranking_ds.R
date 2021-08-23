## Ranking Individuals

# Rank Matrices according to David's Scores


dss <- lapply(m.wl, function(x) compete::ds(x,norm = F)) # each cohort's David's Scores
dss.dfs <- lapply(dss, get_dsdf) # put into dataframe
dss.df.all <- do.call('rbind', Map(cbind, dss.dfs, cohort = LETTERS[1:12])) 

# Plot

davidscores <- ggplot(dss.df.all, aes(x=rank, y=ds, group=cohort)) + 
  geom_line() +
  geom_hline(yintercept = 0, col='red', lty=2) + #for N=10/group, 4.5 is the equivalent of 0 in regular DS.
  ylab("David's Score") +
  xlab("Rank") +
  theme_bw() +
  scale_x_continuous(breaks=1:10)
dev.off()
ggsave("results_figures/davidscores.png", davidscores, width=4,height=4)


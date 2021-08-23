# 1. Make Win-Loss Sociomatrices
df.groupxx <- df.groupx %>% 
  map(~ select(., winner,loser,result)) %>% 
  map(~ as.data.frame(.))

m.wl <- df.groupxx %>% 
  map(get_wl_matrix_new)


# 2. Make Binary Sociomatrices
m.di <- m.wl %>% 
  map(get_di_matrix)


# 3. Calculate modifed h'
m.dv <- lapply(m.wl,devries)
m.dv.p <- lapply(m.dv, function(x) x$`h-modified`) %>% unlist
m.dv.pval <- lapply(m.dv, function(x) x$`p-value`) %>% unlist


# 4. Calculate directional consistency
m.dc <- lapply(m.wl,dc_test)
m.dc.p <- lapply(m.dc, function(x) x$DC) %>% unlist
m.dc.pval <- lapply(m.dc, function(x) x$`DC.pvalue`) %>% unlist


# 5. Steepness
m.st.p <- lapply(m.wl, steepness::getStp) %>% unlist
m.st.pval <- lapply(m.wl, getStp.pval) %>% unlist


# 6. Triangle transitivity
m.tt <- lapply(m.di,ttri_test)
m.tt.p <- lapply(m.tt, function(x) x$ttri) %>% unlist
m.tt.pval <- lapply(m.tt, function(x) x$pval) %>% unlist



# 7. Despotism
m.d <- lapply(m.wl,despotism)
m.d.val <- lapply(m.d, function(x) x[[1]]) %>% unlist
m.d.val


# 8. Gini-coefficients 
gcw <- lapply(m.wl, function(x) ineq::Gini(rowSums(x))) %>% unlist() #GC Wins
gcl <- lapply(m.wl, function(x) ineq::Gini(colSums(x))) %>% unlist() #GC Losses




### Save all results data to results folder

## Store Matrices
matrices <- list(m.wl, m.di)
saveRDS(matrices, "results_statRDS/matrices.RDS")


## Hierarchy Results
data.frame(
  'cohort' = LETTERS[1:12],
  'hvalues' = m.dv.p,
  'dc' = m.dc.p,
  'steepness' = m.st.p,
  'ttri' = m.tt.p,
  'despotism' = m.d.val/100,
  'gini.win' = gcw,
  'gini.lose'= gcl,
  'hvalue.pval' = m.dv.pval,
  'dc.pval' = m.dc.pval,
  'steep.pval' = m.st.pval,
  'ttri.pval' = m.tt.pval
) -> resultsdf

#round
resultsdf[,2:8]  <- round(resultsdf[,2:8],2)
resultsdf[,9:12]  <- round(resultsdf[,9:12],3)
# 
saveRDS(resultsdf, "results_statRDS/resultsdf.RDS")

resultsdf %>% 
  select(-cohort) %>% 
  summarize_all(median_iqr) -> sum

sum[8:length(sum)] <- ""

options(digits = 3)

resultsdfx <- resultsdf %>%
  select(-cohort) %>%
  mutate(hvalue.pval =format(hvalue.pval,nsmall = 3),
         dc.pval =format(dc.pval,nsmall = 3),
         steep.pval =format(steep.pval,nsmall = 3),
         ttri.pval =format(ttri.pval,nsmall = 3))

rbind(resultsdfx %>% 
        mutate_all(as.character),
      sum) %>% 
  mutate(hvalues = glue('{hvalues} ({hvalue.pval})'),         
         dc = glue('{dc} ({dc.pval})'),
         steepness = glue('{steepness} ({steep.pval})'),
         ttri = glue('{ttri} ({ttri.pval})')) -> supp_table


rownames(supp_table) <- c(LETTERS[1:12],"Median (IQR)")
  
write.csv(supp_table[1:7], "results_tables/supp_table.csv", row.names = T)


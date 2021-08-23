


# install.packages('CoDiNA')
# install.packages('wTO')
library(CoDiNA)
library(wTO)

# my data ==================

my_tissue = "Liver"
my_tissue = "Spleen"

# datExpr <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
# 
# rnaseq_sampleid %>% 
#   left_join(all_behavior %>% select(subjectID, domgroup)) %>% 
#   filter(domgroup == "Alpha" & tissue == my_tissue) %>% 
#   .$sampleID -> alphaID
# 
# rnaseq_sampleid %>% 
#   left_join(all_behavior %>% select(subjectID, domgroup)) %>% 
#   filter(domgroup != "Alpha" & tissue == my_tissue) %>% 
#   .$sampleID -> subID
# 
# datExpr[rownames(datExpr) %in% alphaID,] %>% 
#   t %>%  
#   as.data.frame %>% 
#   rownames_to_column('ensgene') %>% 
#   left_join(grcm38 %>% select(ensgene, symbol)) %>% 
#   select(-ensgene) %>% 
#   filter(!is.na(symbol)) %>% 
#   filter(!duplicated(symbol)) %>% 
#   column_to_rownames('symbol')-> alpha_expression
# 
# datExpr[rownames(datExpr) %in% subID,] %>% 
#   t %>%  
#   as.data.frame %>% 
#   rownames_to_column('ensgene') %>% 
#   left_join(grcm38 %>% select(ensgene, symbol)) %>% 
#   select(-ensgene) %>% 
#   filter(!is.na(symbol)) %>% 
#   filter(!duplicated(symbol)) %>% 
#   column_to_rownames('symbol')-> sub_expression



# wTO_out_alpha = wTO.fast(Data = alpha_expression, n = 100)
# saveRDS(wTO_out_alpha,glue("results_RNAseqRDS/wTO_out_alpha_{my_tissue}.RDS"))
wTO_out_alpha <- readRDS(glue("results_RNAseqRDS/wTO_out_alpha_{my_tissue}.RDS"))

# wTO_out_sub = wTO.fast(Data = sub_expression, n = 100)
# saveRDS(wTO_out_sub,glue("results_RNAseqRDS/wTO_out_sub_{my_tissue}.RDS"))
wTO_out_sub <- readRDS(glue("results_RNAseqRDS/wTO_out_sub_{my_tissue}.RDS"))

# wTO_out_alpha %>% 
#   rename(HEAD = Node.1,
#          TAIL = Node.2,
#          WEIGHT = wTO) %>% 
#   select(HEAD, TAIL, WEIGHT) %>% 
#   write.table(.,"network.txt",
#               row.names = F, quote = F, sep = "\t")


# explore ============================================================
# wTO_out_alpha %>% 
#   ggplot(aes(-log10(pval)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Alpha wTO unadjusted p-value"))
# 
# wTO_out_sub %>% 
#   ggplot(aes(-log10(pval)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Sub wTO unadjusted p-value"))
# 
# wTO_out_alpha %>% 
#   ggplot(aes(-log10(pval.adj)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Alpha wTO adjusted p-value"))
# 
# wTO_out_sub %>% 
#   ggplot(aes(-log10(pval.adj)))+
#   geom_histogram()+
#   geom_vline(aes(xintercept = -log10(0.05)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: Sub wTO adjusted p-value"))
# 
# 
# # not so pleasant looking  
# 
# wTO_out_alpha %>% 
#   filter(pval.adj <0.05) %>% 
#   nrow()
# 
# wTO_out_sub %>% 
#   filter(pval.adj <0.05) %>% 
#   nrow()
# 
# 
# wTO_out_sub %>% 
#   arrange(pval)
# 
# # ALPHA = subset(wTO_out_alpha, wTO_out_alpha$pval.adj <0.05, select = c('Node.1', 'Node.2', 'wTO'))
# # SUB = subset(wTO_out_sub, wTO_out_sub$pval.adj <0.05, select = c('Node.1', 'Node.2', 'wTO'))
# 
# 
# wTO_out_alpha %>% 
#   arrange(desc(wTO))


wTO_out_alpha %>% 
  mutate(pval = replace(pval,pval == 0,0.001)) %>% 
  mutate(wTO_logpval = sin(wTO)*-log10(pval)) -> wTO_out_alpha

wTO_out_sub %>% 
  mutate(pval = replace(pval,pval == 0,0.001)) %>% 
  mutate(wTO_logpval = sin(wTO)*-log10(pval)) -> wTO_out_sub


# top_frac(wTO_out_alpha %>% mutate(wTO_2 = abs(wTO)),0.001) %>%
#   arrange(abs(wTO)) %>%
#   head()
# 
# top_frac(wTO_out_sub %>% mutate(wTO_2 = abs(wTO)),0.001) %>%
#   arrange(abs(wTO)) %>%
#   head()


# top_frac(rbind(wTO_out_alpha, wTO_out_sub) %>% mutate(wTO_2 = abs(wTO)),0.01) %>%
#   arrange(abs(wTO)) %>%
#   head()



# 
# top_frac(wTO_out_sub,0.001) %>% 
#   arrange(wTO_logpval) %>% 
#   head
# 
# top_frac(rbind(wTO_out_alpha, wTO_out_sub),0.001) %>% 
#   arrange(wTO_logpval) %>% 
#   head
# 



# top_frac(wTO_out_alpha,0.001) %>%
#   arrange(wTO_logpval) %>% 
#   head
# 
# top_frac(wTO_out_sub,0.001) %>% 
#   arrange(wTO_logpval) %>% 
#   head
# 
# top_frac(rbind(wTO_out_alpha, wTO_out_sub),0.001) %>% 
#   arrange(wTO_logpval) %>% 
#   head
# 

if(my_tissue == "Liver"){my_wTO_cutoff = 0.378} else{my_wTO_cutoff = 0.794}

# my_cutoff = 0.57 # liver
# my_cutoff = 0.5402467 # liver 2 - 0.001 top_frac
# my_cutoff = 2.15 # spleen
# my_cutoff = 2.156244# spleen 2

# wTO_out_alpha %>% 
#   mutate(status = "Alpha") %>% 
#   rbind(wTO_out_sub %>% mutate(status = "Subordinate")) %>% 
#   mutate(pval = replace(pval,pval == 0,0.001)) %>% 
#   mutate(wTO_logpval = sin(wTO)*-log10(pval)) %>% 
#   ggplot(aes(wTO_logpval))+
#   geom_histogram()+ # interesting 
#   geom_vline(xintercept = c(my_cutoff,-(my_cutoff)), color = "red", linetype = "dashed")+
#   theme_bw(base_size = 15)+
#   labs(title = glue("{my_tissue}: sin(wTO)*-log10(pval)"),
#        x = "sin(wTO)*-log10(pval)")+
#   facet_wrap(~status)
# 


# my arbitrary cutoff 
# ALPHA = subset(wTO_out_alpha, abs(wTO_out_alpha$wTO_logpval) > my_cutoff, select = c('Node.1', 'Node.2', 'wTO'))
# SUB = subset(wTO_out_sub, abs(wTO_out_sub$wTO_logpval) > my_cutoff, select = c('Node.1', 'Node.2', 'wTO'))

# my wTO cutoff 
ALPHA = subset(wTO_out_alpha, abs(wTO_out_alpha$wTO) > my_wTO_cutoff, select = c('Node.1', 'Node.2', 'wTO'))
SUB = subset(wTO_out_sub, abs(wTO_out_sub$wTO) > my_wTO_cutoff, select = c('Node.1', 'Node.2', 'wTO'))


# pvalue cutoff
# ALPHA = subset(wTO_out_alpha, wTO_out_alpha$pval <0.05, select = c('Node.1', 'Node.2', 'wTO'))
# SUB = subset(wTO_out_sub, wTO_out_alpha$pval <0.05, select = c('Node.1', 'Node.2', 'wTO'))




DiffNet = MakeDiffNet(Data = list(ALPHA, SUB), Code = c('ALPHA', 'SUB'))

saveRDS(DiffNet, glue("results_RNAseqRDS/{my_tissue}_final_DiffNet.RDS"))





# 
# 
# Graph = plot.CoDiNA(DiffNet, 
#                     # layout = 'layout_components', 
#                     path = glue("FinalVis_max_{my_tissue}_cutoff{my_wTO_cutoff}wTO.html"))
# 
# # str(DiffNet)
# # View(DiffNet)
# 
# # c(DiffNet$Node.1,DiffNet$Node.2) %>% unique() %>% length() -> total_nodes
# nrow(DiffNet) -> total_links
# 
# table(DiffNet$Phi_tilde) %>% as.data.frame()
# 
# # Using the median
# # int_C = quantile(DiffNet$Score_internal, 0.5)
# # ext_C = quantile(DiffNet$Score_Phi, 0.5)
# 
# # # Using the first and the third quantile
# int_C = quantile(DiffNet$Score_internal, 0.25) # figure with 2 is with this option 
# ext_C = quantile(DiffNet$Score_Phi, 0.75)
# 
# Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = ext_C, cutoff.internal = int_C)
# # Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = 0, cutoff.internal = 1)
# 
# nrow(Nodes_Groups) -> total_nodes
# table(Nodes_Groups$Phi_tilde)
# 
# 
# table(Nodes_Groups$Phi_tilde)/nrow(Nodes_Groups) 
# 
# my_tissue
# DiffNet
# round(table(Nodes_Groups$Phi_tilde)/nrow(Nodes_Groups)*100,2)
# round(table(DiffNet$Phi_tilde)/total_links*100 ,2)
# 
# 
# table(DiffNet$Phi_tilde) %>% 
#   as.data.frame() %>% 
#   mutate(perc = Freq/total_links) %>% 
#   mutate(xx = "Links") -> links 
# table(Nodes_Groups$Phi_tilde) %>% 
#   as.data.frame() %>% 
#   mutate(perc = Freq/total_nodes) %>% 
#   mutate(xx = "Nodes") -> nodes 
# 
# rbind(nodes, links) %>% 
#   mutate(xx = factor( xx , levels = c("Nodes", "Links"))) %>% 
#   mutate(Var1 = factor(Var1, levels = rev(c("a","g.ALPHA","g.SUB","U")))) %>% 
#   ggplot(aes(x = xx, y =  perc, fill = Var1))+
#   geom_col(position = "stack",color ='white', alpha = 0.8)+
#   scale_fill_manual(values = rev(c("#a55341", "Purple4","Orange","grey70")),
#                     label = rev(c("Both","Alpha","Subordinate","Undefined")))+
#   theme_minimal(base_size = 8)+
#   labs(x = "",
#        y = "",
#        fill = "Nodes/links specific to:")+
#   theme(panel.grid = element_blank(),
#         axis.text.y = element_text(hjust = 1.5, size = 10),
#         legend.position = "top")+
#   coord_flip() -> linknode_stacked
# 
# linknode_stacked
# 
# png(filename = glue("results_figures/linknode_stacked_{my_tissue}wTO.png"),
#     width = 14, height = 5, units = "cm", res = 600)
# print(linknode_stacked)
# invisible(dev.off())
# 
# 
# 
# 
# 
# 
# 
# 
# Graph = plot(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, 
#             layout = 'layout_components', path = glue("Vis_{my_tissue}_cutoff{my_wTO_cutoff}wTO.html"))
# 
# 
# Graph = plot(DiffNet,cutoff.external = ext_C, cutoff.internal = int_C, 
#              # layout = 'layout_components', 
#              path = glue("Vis_max_{my_tissue}_cutoff{my_wTO_cutoff}wTO.html"))
# 
# 
# 
# Graph = plot.CoDiNA(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, 
#              layout = 'layout_components', path = glue("newcol_Vis_{my_tissue}_cutoff{my_wTO_cutoff}wTO.html"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# DiffNet %>% 
#   filter(Node.1 %in% focus) %>% 
#   filter(Node.2 %in% focus) -> DiffNet2
# 
# 
# 
# 
# Graph = plot(DiffNet2, cutoff.external = ext_C, cutoff.internal = int_C, 
#              layout = 'layout_components', path = glue("Vis_DiffNet2_{my_tissue}_cutoff{my_cutoff}.html"))
# 
# 
# 
# 
# 
# 
# 
# Graph = plot(DiffNet, cutoff.external = 0, cutoff.internal = 1, 
#              layout = 'layout_components', path = glue("Vis_{my_tissue}_cutoff{my_cutoff}.html"))
# # Graph = plot.CoDiNA(DiffNet, cutoff.external = ext_C, cutoff.internal = int_C, 
# #                     layout = 'layout_components', path = glue("Vis_{my_tissue}_cutoff{my_cutoff}2.html"))
# # 
# # Graph
# # 
# # grcm38 %>% 
# #   filter(ensgene == "ENSMUSG00000022877") %>% 
# #   select(symbol, description)
# # 
# # grcm38 %>% 
# #   filter(ensgene == "ENSMUSG00000031765") %>% 
# #   select(symbol, description)
# 
# glimpse(Graph)
# 
# 
# my_top_selection = 100
# Graph$Nodes %>% 
#   select(id,Phi_tilde,Degree_Total) %>% 
#   top_n(my_top_selection) %>% 
#   select(id, Phi_tilde,Degree_Total) -> nodes
# 
# nodes$id -> focus
# focus
# 
# Graph$Edges %>% 
#   filter(id1 %in% focus) %>% 
#   filter(id2 %in% focus) %>% 
#   select(id1, id2,Score)-> edges
# 
# nodes %>% 
#   group_by(Phi_tilde) %>% 
#   count()
# 
# 
# Graph$Nodes %>% 
#   select(id,Phi_tilde,Degree_Total) %>% 
#   top_n(10) %>% 
#   .$id -> labelyes
# 
# 
# 
# net <- igraph::graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE) %>% 
#   igraph::simplify(., remove.multiple = T) %>% 
#   ggnetwork::ggnetwork(.) %>% 
#   mutate(labelyes = ifelse(name %in% labelyes, as.character(name),""))
# 
# 
# 
# 
# ggplot(net,aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey50",alpha=0.2)+
#   geom_nodes(aes(color = Phi_tilde, size = Degree_Total),alpha = 0.9)+
#   geom_nodetext_repel(aes(label = labelyes),
#                 fontface = "bold",
#                 size =1.5) +
#   scale_color_manual(values = c("Purple4","Orange","grey70")) + # for liver
#   # scale_color_manual(values = c("#a55341", "Purple4","Orange","grey70")) + # for spleen
#   theme_blank()+
#   theme(legend.position = "none")  -> my_network_plot
# 
# my_network_plot
# 
# 
# png(filename = glue("results_figures/{my_tissue}_CoDiNA_top{my_top_selection}.png"),
#     width = 10, height = 7, units = "cm", res = 600)
# 
# print(my_network_plot)
# 
# dev.off()
# 
# # just legend - do with Spleen data  =======================================================
# 
# ggplot(net,aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_edges(color = "grey50",alpha=0.2)+
#   geom_nodes(aes(color = Phi_tilde),alpha = 0.9, size = 8)+
#   scale_color_manual(values = c("#a55341", "Purple4","Orange","grey70"),
#                      label = c("Both","Alpha","Subordinate","Undefined"))+
#   theme(legend.position = "top")+
#   labs(color = "Nodes specific to:") -> for_legend
# 
# 
# png(filename = glue("results_figures/legend_CoDiNA.png"),
#     width = 13, height = 2, units = "cm", res = 600)
# 
# leg <- get_legend(for_legend)
# grid.arrange(leg)
# 
# dev.off()
# 

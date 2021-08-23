# Liver 

my_top_selection = 50

Graph <- readRDS("results_RNAseqRDS/CoDiNA_Graph_Liver.RDS")

Graph$Nodes %>% 
  select(-Degree_a,-Degree_b,-Degree_g) %>% 
  top_n(my_top_selection) %>%
  select(id, Phi_tilde,Degree_Total,Phi) %>% 
  .$id -> focus

DiffNet <- readRDS("results_RNAseqRDS/Liver_final_DiffNet.RDS")


Graph = plot.CoDiNA.og(DiffNet, 
                       cutoff.external =  quantile(DiffNet$Score_Phi, 0.75),
                       cutoff.internal = quantile(DiffNet$Score_internal, 0.25), 
                       layout = 'layout_nicely', 
                       path = glue("FinalVis_OG_Liver.html"))

# Graph$Nodes %>% nrow() -> my_top_selection

Graph$Nodes %>% 
  select(-Degree_a,-Degree_b,-Degree_g) %>% 
  top_n(my_top_selection) %>%
  select(id, Phi_tilde,Degree_Total,Phi) %>% 
  .$id -> focus

focus
# Graph$Nodes %>% 
#   select(id,Phi_tilde,Degree_Total) %>% 
#   group_by(Phi_tilde) %>% 
#   top_n(15) %>% 
#   filter(Phi_tilde != "U") %>% 
#   filter(id != "2200002D01Rik") %>% 
#   .$id -> labelyes


newGraph = plot.CoDiNA(DiffNet %>% 
                         filter(Node.1 %in% focus) %>%
                         filter(Node.2 %in% focus), 
                       cutoff.external =  quantile(DiffNet$Score_Phi, 0.75),
                       cutoff.internal = quantile(DiffNet$Score_internal, 0.25), 
                       layout = 'layout_nicely', 
                       cluster = F,
                       path = glue("filter_FinalVis_Liver_topgene{my_top_selection}wTO.html"))

newGraph


Graphx = plot.CoDiNA(DiffNet, 
                    cutoff.external =  quantile(DiffNet$Score_Phi, 0.75),
                    cutoff.internal = quantile(DiffNet$Score_internal, 0.25), 
                    layout = 'layout_nicely', 
                    path = glue("FinalVis_Liver.html"))


Graphx




# Stacked bar

# c(DiffNet$Node.1,DiffNet$Node.2) %>% unique() %>% length() -> total_nodes
nrow(DiffNet) -> total_links

table(DiffNet$Phi_tilde) %>% as.data.frame()

# Using the median
# int_C = quantile(DiffNet$Score_internal, 0.5)
# ext_C = quantile(DiffNet$Score_Phi, 0.5)

# # Using the first and the third quantile
int_C = quantile(DiffNet$Score_internal, 0.25) # figure with 2 is with this option
ext_C = quantile(DiffNet$Score_Phi, 0.75)

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = ext_C, cutoff.internal = int_C)
# Nodes_Groups = ClusterNodes(DiffNet = DiffNet, cutoff.external = 0, cutoff.internal = 1)

nrow(Nodes_Groups) -> total_nodes
table(Nodes_Groups$Phi_tilde)


table(Nodes_Groups$Phi_tilde)/nrow(Nodes_Groups)

my_tissue = "Liver"
DiffNet
round(table(Nodes_Groups$Phi_tilde)/nrow(Nodes_Groups)*100,2)
round(table(DiffNet$Phi_tilde)/total_links*100 ,2)


table(DiffNet$Phi_tilde) %>%
  as.data.frame() %>%
  mutate(perc = Freq/total_links) %>%
  mutate(xx = "Links") -> links
table(Nodes_Groups$Phi_tilde) %>%
  as.data.frame() %>%
  mutate(perc = Freq/total_nodes) %>%
  mutate(xx = "Nodes") -> nodes

rbind(nodes, links) %>%
  mutate(xx = factor( xx , levels = c("Nodes", "Links"))) %>%
  mutate(Var1 = factor(Var1, levels = rev(c("a","g.ALPHA","g.SUB","U")))) %>%
  ggplot(aes(x = xx, y =  perc, fill = Var1))+
  geom_col(position = "stack",color ='white', alpha = 0.8)+
  scale_fill_manual(values = rev(c("limegreen", "Purple4","Orange","grey70")),
                    label = rev(c("Both","Alpha","Subordinate","Undefined")))+
  theme_minimal(base_size = 8)+
  labs(x = "",
       y = "",
       fill = "Nodes/links specific to:")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(hjust = 1.5, size = 10),
        legend.position = "top")+
  coord_flip() -> linknode_stacked

linknode_stacked

png(filename = glue("results_figures/linknode_stacked_{my_tissue}wTOnewcol.png"),
    width = 14, height = 5, units = "cm", res = 600)
print(linknode_stacked)
invisible(dev.off())


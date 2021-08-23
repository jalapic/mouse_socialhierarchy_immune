# https://github.com/cran/CoDiNA/blob/master/R/plot2.R

DiffNet %>% 
  filter(Node.1 %in% focus) %>%
  filter(Node.2 %in% focus) -> x

x <- DiffNet
cutoff.external = ext_C
cutoff.internal = int_C
cutoff.ratio = 1
# layout = 'layout_components'
layout = "layout_nicely"
smooth.edges = TRUE
path = NULL
MakeGroups = FALSE
Cluster = T
legend = TRUE
manipulation = FALSE
sort.by.Phi = T


plot.CoDiNA = function(x, cutoff.external = 0, cutoff.internal = 1, cutoff.ratio = 1,
                       layout = NULL, smooth.edges = TRUE,
                       path = NULL, MakeGroups = FALSE,
                       Cluster = FALSE, legend = TRUE,
                       manipulation = FALSE,
                       sort.by.Phi = FALSE , ...)
{
  `%ni%` <- Negate(`%in%`)
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  clean = subset(x, x$Score_Phi_tilde > cutoff.external &
                   x$Score_internal < cutoff.internal &
                   x$Score_ratio > cutoff.ratio)
  
  Vars = c('Node.1', 'Node.2', 'Score_Phi_tilde',
           'Score_internal', 'Phi', 'Phi_tilde', 'Score_ratio')
  # message(Vars)
  if (any(Vars %ni% names(x))) {
    stop("x input is not complete.")
  }
  Score = by =  'Score_Phi_tilde'
  
  input_vis = data.frame(subset(clean,!is.na(clean$Node.1) ,select = Vars))
  
  
  if (is.numeric(cutoff.external) == FALSE | is.numeric(cutoff.internal) == FALSE) {
    stop("cutoff value must be numeric.")
  }
  if (Cluster %ni% c(TRUE, FALSE)) {
    stop("Cluster must be TRUE or FALSE.")
  }
  if (smooth.edges %ni% c(TRUE, FALSE)) {
    stop("smooth.edges must be TRUE or FALSE.")
  }
  if(cutoff.external<0.01){
    input_vis = droplevels(subset(input_vis, abs(input_vis$Score_Phi_tilde) > 0.01))
  }
  if (nrow(input_vis) <= 2) {
    stop("Not enough nodes on your network. Choose a another cutoff.")
  }
  if (smooth.edges == TRUE) {
    smooth.edges = "enabled"
  }
  
  Nodes_Phi = ClusterNodes(DiffNet = x, cutoff.external = cutoff.external, cutoff.internal = cutoff.internal)
  
  input_vis = input_vis[!is.na(input_vis$Score_Phi_tilde), ]
  input_vis = input_vis[!is.na(input_vis$Score_internal), ]
  
  ####
  #### getting colors for Phi groups
  ####
  input_vis$GROUP_FULL = as.factor(input_vis$Phi_tilde)
  
  
  colsC = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'a'))))
  colsS = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'g'))))
  colsD = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'b'))))
  
  
  
  if(nrow(colsD) == 0){
    colsD = data.frame(GROUP_FULL = NULL,
                       Phi = NULL, color = NULL, shape = NULL)
  }
  if(nrow(colsS) == 0){
    colsD = data.frame(GROUP_FULL = NULL, Phi = NULL,color = NULL, shape = NULL)
  }
  if(nrow(colsC) == 0){
    colsC = data.frame(GROUP_FULL = NULL, Phi = NULL,color = NULL, shape = NULL)
  }
  if(nrow(colsC) > 0){
    colsC$color = colorRampPalette(c('limegreen', 'olivedrab1'))(nrow(colsC))
    colsC$shape = 'triangle'
    colsC$Phi = 'a'
  }
  if(nrow(colsS) > 0){
    colsS$color = colorRampPalette(c('#6c1bb5', 'orange'))(nrow(colsS))
    colsS$shape = 'star'
    colsS$Phi = 'g'
  }
  if(nrow(colsD) > 0){
    colsD$color = colorRampPalette(c('salmon','violetred3'))(nrow(colsD))
    colsD$shape = 'square'
    colsD$Phi = 'b'
  }
  
  colsI = data.frame(GROUP_FULL = 'U', color = 'lightgrey',
                     shape = 'diamond', Phi = 'U')
  
  colormap = rbind(colsC, colsD, colsS, colsI)
  input_vis = suppressMessages(plyr::join(input_vis, colormap))
  
  input_vis = suppressMessages(plyr::arrange(input_vis, input_vis$Node.1, input_vis$Node.2))
  input_vis = droplevels(input_vis)
  nodes <- data.frame(id = sort(unique(c(as.character(input_vis$Node.1),
                                         as.character(input_vis$Node.2)))))
  gg = igraph::graph_from_data_frame(data.frame(input_vis$Node.1, input_vis$Node.2, weights = input_vis$Score_Phi_tilde), directed = FALSE)
  
  DEGREE = as.data.frame(igraph::degree(gg))
  if(MakeGroups == FALSE){
    group = 1
  }
  # if (MakeGroups == 'infomap'){
  #   
  #   group = igraph::cluster_infomap(gg)$membership
  # }
  # else if (MakeGroups == 'walktrap'){
  #   group = igraph::cluster_walktrap(gg)$membership
  # }
  # else if (MakeGroups == 'leading_eigen'){
  #   group = igraph::cluster_leading_eigen(gg)$membership
  # }
  # else if (MakeGroups == 'louvain'){
  #   group = igraph::cluster_louvain(gg)$membership
  # }
  # else if (MakeGroups == 'label_prop'){
  #   group = igraph::cluster_label_prop(gg)$membership
  # }
  # else if (MakeGroups == 'fast_greedy'){
  #   group = igraph::cluster_fast_greedy(gg)$membership
  # }
  # else if (MakeGroups == 'optimal'){
  #   group = igraph::cluster_optimal(gg)$membership
  # }
  # else if (MakeGroups == 'spinglass'){
  #   group = igraph::cluster_spinglass(gg)$membership
  # }
  # else if (MakeGroups == 'edge.betweenness'){
  #   group = igraph::edge.betweenness.community(gg)$membership
  # }
  
  
  nodes = suppressMessages(plyr::join(nodes, data.frame(id = igraph::V(gg)$name, cluster = group)))
  
  igraph::E(gg)$weight = abs(input_vis$Score_Phi_tilde)
  names(DEGREE) = "degree"
  DEGREE$id = row.names(DEGREE)
  nodes = suppressMessages(plyr::join(nodes, DEGREE))
  nodes$value = (nodes$degree - min(nodes$degree))/(max(nodes$degree) -
                                                      min(nodes$degree))
  nodes$value = nodes$value^2 + 5
  igraph::V(gg)$size = nodes$value
  igraph::V(gg)$color = nodes$cluster
  
  nodes$size = nodes$value
  names(Nodes_Phi) [1]= 'id'
  nodes = suppressWarnings(suppressMessages( plyr::join(nodes, Nodes_Phi)))
  
  Nodes1 = data.frame(id = clean$Node.1, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)
  Nodes2 = data.frame(id = clean$Node.2, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)
  Nodes = rbind(Nodes1, Nodes2)
  
  
  Map2 =  reshape2::dcast(Nodes, id~Phi, fun.aggregate = length,
                          value.var = 'Phi')
  rm(Nodes1)
  rm(Nodes2)
  rm(Nodes)
  
  
  
  names(Map2)[1] = 'id'
  if(is.null(Map2$a)){
    Map2$a = 0
  }
  if(is.null(Map2$g)){
    Map2$g = 0
  }
  if(is.null(Map2$b)){
    Map2$b = 0
  }
  
  nodes = suppressMessages(plyr::join(nodes, Map2))
  
  nodes$label = nodes$id
  # nodes$label = ifelse(nodes$id %in% labelyes, nodes$id, "")
  
  
  nodes$title = paste0("<p> Node ID: ", nodes$label,
                       "<br>Degree: ", nodes$degree,
                       "<br>Degree a: ", nodes$a,
                       "<br>Degree b: ", nodes$b,
                       "<br>Degree g: ", nodes$g,
                       '<br>Phi:', nodes$Phi,
                       '<br>Phi_tilde:', nodes$Phi_tilde,
                       "</p>")
  if(sort.by.Phi == TRUE){
    nodes$id = paste(nodes$groupPhi, nodes$group, nodes$id, sep ='_')
    
  }
  
  nodes$GROUP_FULL <- nodes$Phi_tilde <- as.character(nodes$Phi_tilde)
  
  # nodes$color = ifelse(nodes$Phi == 'a',  'green', 'gray30')
  # nodes$color = ifelse(nodes$Phi == 'b',  'red', nodes$color)
  # nodes$color = ifelse(nodes$Phi == 'g',  'blue', nodes$color)
  nodes$color = ifelse(nodes$Phi_tilde == 'a',  'green', 'lightgrey')
  nodes$color = ifelse(nodes$Phi_tilde == 'b.ALPHA',  'red', nodes$color)
  nodes$color = ifelse(nodes$Phi_tilde == 'b.SUB',  'red', nodes$color)
  nodes$color = ifelse(nodes$Phi_tilde == 'g.ALPHA',  '#6c1bb5', nodes$color)
  nodes$color = ifelse(nodes$Phi_tilde == 'g.SUB',  'orange', nodes$color)
  table(nodes$color)
    nodes$frame.color = nodes$color
  
  nodes = nodes[order(nodes$id),]
  node1_ID =  data.frame(Node.1 = nodes$label, from = nodes$id)
  node2_ID = data.frame(Node.2 = nodes$label, to = nodes$id)
  # input2 = data.frame(Node.1 = input_vis$Node.1, Node.2= input_vis$Node.2)
  edges_ID = suppressMessages(plyr::join_all(list(input_vis, node1_ID, node2_ID) ))
  edges_ID$label = apply(edges_ID[,1:2], 1, paste, collapse = '<->')
  edges_ID$L1 = edges_ID[,1]
  edges_ID$L2 = edges_ID[,2]
  edges <- data.frame(from = edges_ID$from, to = edges_ID$to,
                      Label = edges_ID$label, L1 = edges_ID$L1, L2 = edges_ID$L2,
                      group = edges_ID$Phi_tilde,
                      Score = edges_ID$Score_Phi_tilde,
                      Phi = edges_ID$Phi, color = edges_ID$color)
  wto = abs(edges_ID$Score_ratio)
  edges$width = 3*abs((wto - min(wto))/(max(wto) -
                                          min(wto)))^4
  nodes = droplevels(nodes)
  edges = droplevels(edges)
  ledges <- data.frame(color = colormap$color,
                       label = colormap$GROUP_FULL,
                       arrows = '')
  ledges2 <- data.frame(color = c('green', 'red', 'red', '#6c1bb5', 'orange', 'lightgrey'),
                        label = c('a', 'b.ALPHA', 'b.SUB', 'g.ALPHA', 'g.SUB', 'U'),
                        arrows = c("", "", '', '','',''))
  ledges2 = rbind(ledges, ledges2)
  edges$title = paste0("<p>Edge: ", edges$Label,
                       "<br>Score: ", round(edges$Score, 2),
                       '<br>Group:', edges$group,
                       '<br>Phi:', edges$Phi,"</p>")
  
  NODESIZE = length(nodes$id)
  EDGESIZE = length(edges$Label)
  main = paste('Network contains:', NODESIZE, 'nodes and', EDGESIZE, 'edges.')
  
  nodes$color.border= nodes$color
  
  
  network <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visInteraction(navigationButtons = FALSE, hover = TRUE, multiselect = FALSE) %>%
    # visNetwork::visEdges(smooth = smooth.edges )%>%
    # visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
    #                                                degree = 1, hover = FALSE),
    #                        nodesIdSelection = list(enabled = TRUE, useLabels =TRUE,
    #                                                style = "width: 200px; height: 26px;\\n   background: #f8f8f8;\\n   color: darkblue;\\n   border:none;\\n   outline:none;"),
    #                        manipulation = F,
    #                        # selectedBy = list(variable = 'cluster', multiple = FALSE),
    #                        selectedBy = list(variable = 'Phi_tilde', multiple = FALSE),
    #                        # collapse = list(enabled = TRUE, clusterOptions =list(Phi = nodes$groupPhi),resetHighlight = TRUE)) %>%
    # ) %>%
    visNetwork::visPhysics(enabled = F)  %>%
    visNetwork::visIgraphLayout(layout = layout) #%>% 
  # visNetwork::visExport(type = "png",
  #                       name = "networkpng",
  #                       float = "right",
  #                       label = "Save png",
  #                       background = "transparent",
  #                       style= "")
  # visNetwork::visSave(network, file = path)
  # message(path)

  return(network)
}

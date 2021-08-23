### Hub Score
funhub<-function(g){
  hg<-hub.score(g)
  return(hg$vector)
}



get_indiv_network_measures <- function(df){
  d1<-get_wl_matrix_new(df %>% select(winner, loser, result))
  
  d2pa<-get_di_matrix(d1, type = "pa")
  d2wl<-get_di_matrix(d1, type = "wl")
  
  gdfpab<-graph.adjacency(d2pa)
  
  #outdegree
  inet.outdeg<-as.data.frame(igraph::degree(gdfpab, mode=c("out"))) 
  colnames(inet.outdeg)<-c("inet.outdeg")
  
  #indegree
  inet.indeg<-as.data.frame(igraph::degree(gdfpab, mode=c("in"))) 
  colnames(inet.indeg)<-c("inet.indeg")
  
  #out-closeness
  inet.outclose<-as.data.frame(igraph::closeness(gdfpab, mode=c("out"), normalized=TRUE))
  colnames(inet.outclose)<-c("inet.outclose")
  
  #incloseness
  inet.inclose<-as.data.frame(igraph::closeness(gdfpab, mode=c("in"), normalized=TRUE))
  colnames(inet.inclose)<-c("inet.inclose")
  
  #betweenness centrality 
  inet.betw<-as.data.frame(igraph::betweenness(gdfpab,  normalized=TRUE))
  colnames(inet.betw)<-c("inet.betw")
  
  ### For power analyses, use win/loss binarized data
  
  gagwl<-graph.adjacency(d2wl)
  
  gamma<- -0.09 #exponent value used in Bonacich's Power Calculations
  
  inet.bonpow<-as.data.frame(igraph::bonpow(gagwl,exponent=gamma, rescale=TRUE))
  colnames(inet.bonpow)<-c("inet.bonpow1")
  
  inet.hub<-as.data.frame(funhub(gagwl))
  colnames(inet.hub)<-c("inet.hub")
  
  
  network_result<-cbind(mouseID=row.names(inet.outdeg),
                        inet.outdeg,
                        inet.indeg,
                        inet.outclose,
                        inet.inclose,
                        inet.betw,
                        inet.bonpow,
                        inet.hub)
  
  return(network_result)}





## to get HWIG 
agg_interaction_freq <- function(df){
  
  x <- df %>% 
    select(winner,loser) %>% 
    get_wl_matrix()
  
  tx <- x + t(x)  
  
  txx <- tx %>% 
    as.data.frame() %>% 
    mutate(mouseA = as.character(winner),
           mouseB = as.character(loser)) %>% 
    select(mouseA,mouseB,Freq)
  
  return(txx)
}





###HWI - to ultimately calculate HWIG (Godde et al. 2013 Animal Behavior)
hwi_fun <- function(df2, grankA=NULL, grankB=NULL){
  
  df2x <- df2 %>% filter(grank==grankA|grank==grankB) %>% select(realtime, side, grank)
  df2z<-data.table(df2x)
  df2Z <- data.table::dcast(df2z, realtime ~ side, 
                            value.var='grank',
                            fill=NA,  
                            fun.aggregate=list)
  
  ddf <- data.frame(time = df2Z[,1])
  ddf$leftx <- apply(df2Z[,2], 1, function(x) grepl(",", x) )
  ddf$rightx <- apply(df2Z[,3], 1, function(x) grepl(",", x) )
  ddf$left1 <-  apply(df2Z[,2], 1, function(x) grepl(paste0("\\b",grankA,"\\b"), x) )
  ddf$right1 <-  apply(df2Z[,3], 1, function(x) grepl(paste0("\\b",grankA,"\\b"), x) )
  ddf$left2 <-  apply(df2Z[,2], 1, function(x) grepl(paste0("\\b",grankB,"\\b"), x) )
  ddf$right2 <-  apply(df2Z[,3], 1, function(x) grepl(paste0("\\b",grankB,"\\b"), x) )
  ddf$sum1 <- ddf$left1 + ddf$right1
  ddf$sum2 <- ddf$left2 + ddf$right2
  
  #rules
  ddf$value <- ifelse(ddf$leftx==T|ddf$rightx==T, "x",
                      ifelse( (ddf$left1==T & ddf$right2==T) |   (ddf$left2==T & ddf$right1==T)  , "yAB",
                              ifelse(ddf$sum1>0 & ddf$sum2 ==0, "yA",
                                     ifelse(ddf$sum2>0 & ddf$sum1==0, "yB", NA))))
  
  #get values
  x <- sum(ddf$value=="x")
  yAB <- sum(ddf$value=="yAB")
  yA <- sum(ddf$value=="yA")
  yB <- sum(ddf$value=="yB")
  
  hwi <- x / (x + yAB + 0.5*yA + 0.5*yB)
  hwi
  
  return(HWI = hwi)
}

hwi_fun_all<- function(df2){
  temp=data.frame(grankA=NA,
                  grankB=NA,
                  HWI=NA)
  for(i in 1:12){
    for(j in 1:12){
      if (i==j){ a=data.frame(grankA=i,grankB=j,HWI=NA)
      temp<-bind_rows(temp,a)}
      else if (i>j){  a=data.frame(grankA=i,grankB=j,HWI=NA)
      temp<-bind_rows(temp,a)}
      else if (i<j){b=data.frame(grankA=i,grankB=j,HWI=hwi_fun(df2,grankA=i,grankB=j))
      temp<-bind_rows(temp,b) 
      }
    }
  }
  result<-temp %>% filter(!is.na(HWI))
  return(result)
}

hwig_fun<-function(hwi_df,i,j){
  temp_i<-hwi_df %>% filter(grankA==i|grankB==i) %>% summarise(sum_hwi=sum(HWI)) 
  temp_j<-hwi_df %>% filter(grankA==j|grankB==j) %>% summarise(sum_hwi=sum(HWI)) 
  temp_ij<-hwi_df %>% filter(grankA==i&grankB==j) %>% select(HWI)
  sum_hwi=sum(hwi_df$HWI)
  
  hwig<-temp_ij[1,1]*sum_hwi/(temp_i[1,1]*temp_j[1,1])
  return(hwig)
}

hwig_fun_all<- function(hwi_df){
  temp=data.frame(grankA=NA,
                  grankB=NA,
                  HWIG=NA)
  for(i in 1:12){
    for(j in 1:12){
      if (i==j){ a=data.frame(grankA=i,grankB=j,HWIG=NA)
      temp<-bind_rows(temp,a)}
      else if (i>j){  a=data.frame(grankA=i,grankB=j,HWIG=NA)
      temp<-bind_rows(temp,a)}
      else if (i<j){b=data.frame(grankA=i,grankB=j,HWIG=hwig_fun(hwi_df,i,j))
      temp<-bind_rows(temp,b) 
      }
    }
  }
  result<-temp %>% filter(!is.na(HWIG))
  return(result)
}

### calculate gregariousness

greg<-function(df){
  temp_i<-list()
  for(i in 1:max(df$grankB)){
    temp_i[[i]]<-df %>% filter(grankA==i|grankB==i) %>% 
      summarise(sum_HWI=sum(HWI)) %>% 
      mutate(glicko_rank=i) %>% 
      as.data.frame()
  }
  temp<-rbindlist(temp_i) %>% as.data.frame()
  return(temp)
}

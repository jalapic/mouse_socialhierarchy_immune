

df.behav_freq
df.behav_perc
df.sub_behav_freq
df.sub_behav_perc

# All behavior, frequency
df <- df.behav_freq %>%
  select(-mouseID, -cohort) %>% 
  column_to_rownames(var = "subjectID") 
# 
# dfx <- df.behav_freq %>%
#   select(-mouseID) %>% 
#   mutate(cohort = as.numeric(as.factor(cohort))) %>% 
#   column_to_rownames(var = "subjectID") 
# 
# 
# 
# # WHAT if I do behavior frequency per observation hour?
# hour
# 
# 
dfxx <- df.behav_freq %>%
  left_join(hour) %>%
  left_join(all_behavior %>% select(subjectID, inet.betw)) %>% 
  select(-mouseID, -cohort) %>%
  mutate_if(is.numeric,~.x/totalobs) %>%
  select(-totalobs) %>%
  column_to_rownames(var = "subjectID")
# 
# 
# # basically same either normalize by hours or not 
# 


dfxxx <- df.behav_freq %>%
  select(-mouseID, -cohort) %>% 
  left_join(all_behavior %>% select(subjectID, inet.betw)) %>% 
  column_to_rownames(var = "subjectID") 


# dfxx -> my_data
# dfx -> my_data

df -> my_data
fviz_nbclust(my_data, kmeans, method = "gap_stat")
fviz_nbclust(my_data, kmeans, method = "wss")
fviz_nbclust(my_data, kmeans, method = "silhouette")

set.seed(33)
km.res3 <- kmeans(my_data, 3, nstart = 25)
set.seed(33)
km.res4 <- kmeans(my_data, 4, nstart = 25)



# Visualize
fviz_cluster(km.res3, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

fviz_cluster(km.res4, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

km.res3$cluster
km.res4$cluster
set.seed(33)
km.res3x <- kmeans(dfxxx, 3, nstart = 25)

fviz_cluster(km.res3x, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())


set.seed(33)
cluster <- data.frame(subjectID = names(km.res3$cluster),
                      km3_cluster = km.res3$cluster,
                      km4_cluster = km.res4$cluster,
                      km3_cluster_net = km.res3x$cluster) %>% 
  mutate(km3_cluster = paste("Cluster",km3_cluster),
         km4_cluster = paste("Cluster",km4_cluster),
         km3_cluster_net = paste("Cluster",km3_cluster_net))
  



cluster_exp
colnames(cluster_exp)

cluster_exp %>% 
  left_join(all_behavior %>% select(subjectID,cohort)) ->df

library(lme4)
library(lmerTest)

summary(lmer(lunge ~ km3_cluster+(1|cohort), data = df))

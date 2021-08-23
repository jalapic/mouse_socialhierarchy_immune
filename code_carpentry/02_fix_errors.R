#### Fixing Errors in Raw Data Files

c105to106 <- cohort105 %>% filter(Error=="c106")
c106to105 <- cohort106 %>% filter(Error=="c105")
c107to108 <- cohort107 %>% filter(Error=="c108")
c108to107 <- cohort108 %>% filter(Error=="c107")
c109to110 <- cohort109 %>% filter(Error=="c110")
c110to109 <- cohort110 %>% filter(Error=="c109")
c111to112 <- cohort111 %>% filter(Error=="c112")
c112to111 <- cohort112 %>% filter(Error=="c111")
c113to114 <- cohort113 %>% filter(Error=="c114")
c114to113 <- cohort114 %>% filter(Error=="c113")
c115to116 <- cohort115 %>% filter(Error=="c116")
c116to115 <- cohort116 %>% filter(Error=="c115")


cohort105 <- cohort105 %>% filter(is.na(Error)) %>% rbind(.,c106to105) %>% arrange(Timestamp)
cohort106 <- cohort106 %>% filter(is.na(Error)) %>% rbind(.,c105to106) %>% arrange(Timestamp)
cohort107 <- cohort107 %>% filter(is.na(Error)) %>% rbind(.,c108to107) %>% arrange(Timestamp)
cohort108 <- cohort108 %>% filter(is.na(Error)) %>% rbind(.,c107to108) %>% arrange(Timestamp)
cohort109 <- cohort109 %>% filter(is.na(Error)) %>% rbind(.,c110to109) %>% arrange(Timestamp)
cohort110 <- cohort110 %>% filter(is.na(Error)) %>% rbind(.,c109to110) %>% arrange(Timestamp)
cohort111 <- cohort111 %>% filter(is.na(Error)) %>% rbind(.,c112to111) %>% arrange(Timestamp)
cohort112 <- cohort112 %>% filter(is.na(Error)) %>% rbind(.,c111to112) %>% arrange(Timestamp)
cohort113 <- cohort113 %>% filter(is.na(Error)) %>% rbind(.,c114to113) %>% arrange(Timestamp)
cohort114 <- cohort114 %>% filter(is.na(Error)) %>% rbind(.,c113to114) %>% arrange(Timestamp)
cohort115 <- cohort115 %>% filter(is.na(Error)) %>% rbind(.,c116to115) %>% arrange(Timestamp)
cohort116 <- cohort116 %>% filter(is.na(Error)) %>% rbind(.,c115to116) %>% arrange(Timestamp)


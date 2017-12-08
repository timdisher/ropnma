source("./analyses/final/rop_import.R")


#Start by looking at the study level data
summary(rop_data_study)

rop_data_study %>% keep(is.factor) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1)) + facet_wrap(~key, scales = "free")


rop_data_study %>% keep(is.numeric) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1)) + facet_wrap(~key, scales = "free")


#Arm level data
rop_data_arm %>% keep(is.factor) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) + facet_wrap(~key, scales = "free")


rop_data_arm %>% keep(is.numeric) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0)) + facet_wrap(~key, scales = "free")


#Generate missing covariate values
library(mice)

#Select relevant variables to keep 
keep = c("method","speculum","scleral_dep","exam_number","avg_ga","avg_pma","avg_bw")
#Filter out repeated studies
filt = c("Zeraati 2015a","Zeraati 2015b","Zeraati 2015c")

rop_miss = rop_data_study  %>% filter(!studlab %in% (filt)) %>% select(studlab,keep)

orig = rop_miss %>% drop_na() %>% select(-c(studlab:exam_number)) %>% mutate(id = 1:15)
test = ampute(orig %>% select(-c(id)), prop = 0.30)

test_imp = mice(test$amp)

summary(orig)
summary(complete(test_imp))

t = complete(test_imp,action = "long", include = FALSE) %>% rename(id = .id) %>% gather(var,value,avg_ga:avg_bw) %>% arrange(.imp)

dropped = test$amp %>% mutate(id = 1:15) %>% gather(var,dr_val,avg_ga:avg_bw)
orig  = orig %>% gather(var,or_val,avg_ga:avg_bw)

graph = left_join(dropped,orig, by = c("id","var")) %>% mutate(dropped = ifelse(is.na(dr_val),"yes","no"),id = as.factor(id)) %>% select(-dr_val)

graph = t %>% left_join(graph, by = c("id","var")) %>% mutate(value = ifelse(dropped == "no",NA,value))

base = graph %>% filter(.imp == 1)

base %>% ggplot(aes(x = id,y = or_val)) + geom_point() + facet_wrap(~var, scales = "free_y") + geom_point(aes(y = value, x = id),colour = "red", data = graph)

#Conclusion here is that imputation is unlikely to retrieve accurate values in this dataset and the risk is too high for a single study to bias an already unstable
# analysis

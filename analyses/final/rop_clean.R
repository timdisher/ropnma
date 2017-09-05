source("./analyses/final/rop_import.R")


#Start by looking at the study level data
summary(rop_data_study)

rop_data_study %>% keep(is.factor) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1)) + facet_wrap(~key, scales = "free")


rop_data_study %>% keep(is.numeric) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1)) + facet_wrap(~key, scales = "free")


#Arm level data
rop_data_arm%>% keep(is.factor) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) + facet_wrap(~key, scales = "free")


rop_data_arm %>% keep(is.numeric) %>% gather() %>% ggplot(aes(value)) + 
  geom_bar()  + theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0)) + facet_wrap(~key, scales = "free")

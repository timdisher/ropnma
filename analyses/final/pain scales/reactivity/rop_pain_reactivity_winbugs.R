#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# install.packages('coda')
# install.packages('R2WinBUGS')
# install.packages('netmeta')
# install.packages('reshape2')

library(coda)
library(R2WinBUGS)
library(tidyverse)
library(netmeta)
library(stargazer)
library(reshape2)
library(forcats)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
source("./analyses/final/rop_explore_pipp.R")

pa_reac$trt_group = fct_infreq(pa_reac$trt_group) #reorders factor from most to least common

pa_reac = droplevels(pa_reac)
#Convert to long format and ensure treatment order is maintained

(pa_reac_input = pa_reac %>% select(studlab,trt_group,mean,std_dev,sample_size) %>% 
  rename(y = mean,
         n = sample_size) %>%  
    arrange(studlab,trt_group) %>% 
    group_by(studlab) %>%
    mutate(arm = row_number(),
           t = as.numeric(trt_group),
           na = n())
          )

treatments = pa_reac_input %>% ungroup() %>% select(trt_group,t) %>% distinct() %>% arrange(t) %>%
  rename(description = trt_group)

(pa_reac_wide= pa_reac_input %>% select(-trt_group) %>%
  recast(studlab ~ variable + arm, id.var = c("studlab","arm")) %>% select(studlab:na_1) %>% rename(na = na_1)
  )

(pa_reac_wide = pa_reac_wide %>% mutate(y_2 = y_2 - y_1,
                        y_3 = y_3 - y_1,
                        y_4 = y_4 - y_1,
                        se_2 = sqrt((std_dev_1^2*(n_1-1) + std_dev_2^2*(n_2-1))/(n_1+n_2-2))*sqrt(1/n_1+1/n_2),
                        se_3 = sqrt((std_dev_1^2*(n_1-1) + std_dev_3^2*(n_3-1))/(n_1+n_3-2))*sqrt(1/n_1+1/n_3),
                        se_4 = sqrt((std_dev_1^2*(n_1-1) + std_dev_4^2*(n_4-1))/(n_1+n_4-2))*sqrt(1/n_1+1/n_4),
                        V = (std_dev_1/sqrt(n_1))^2) %>% arrange(na))

(pa_reac_wb = pa_reac_wide %>% select(t_1:t_4,y_2:y_4,se_2:se_4,V,na))


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



fe.model = nma_cont_fe(pa_reac_wb,treatments)



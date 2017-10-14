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
library(scales)
library(forestplot)
library(gemtc)
library(reshape2)
library(R2WinBUGS)

source("./analyses/final/spo2/rop_explore_spo2.R")

source("./functions/nma_cont.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
source("./functions/nma_utility_functions.R")
# #========================================================================================
# 
# 
# Outcome: Sp02 Reactivity-----
# 
# #========================================================================================




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
os_reac_data = NULL

os_reac_data$pa$gemtc = prep_gem(os_reac$data)

os_reac_data$pa$gemtc$data = os_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(os_reac$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))



os_reac_data$pa$network = mtc.network(data.re =os_reac_data$pa$gemtc$data[,1:4],
                                       studies =os_reac_data$pa$gemtc$data[,c(1,5:10)])



os_reac_data$pa$results = mtc.model(os_reac_data$pa$network, type = "consistency",
                                     linearModel = "random",likelihood = "normal",
                                     link = "identity")
os_reac_data$pa$results = mtc.run(os_reac_data$pa$results)

summary(os_reac_data$pa$results)


forest(relative.effect.table(os_reac_data$pa$results),"drops")


#====
# Informative priors = 3 points large diff = 1/3^2 = 0.11

os_reac_data$sa1$results = mtc.model(os_reac_data$pa$network, type = "consistency",
                                    linearModel = "random",likelihood = "normal",
                                    link = "identity",
                                    hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.11))
os_reac_data$sa1$results = mtc.run(os_reac_data$sa1$results)

summary(os_reac_data$sa1$results)


forest(relative.effect.table(os_reac_data$sa1$results),"drops")


save(os_reac, file = "./cache/spo2_reac.rda")


# #========================================================================================
# 
# 
# Outcome: Sp02 Recovery-----
# 
# #========================================================================================

##==== Single study

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#==== Single study only and treatments are equivalent
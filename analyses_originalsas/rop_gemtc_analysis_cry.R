#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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
library(coda)


source("./analyses/final/cry time/rop_explore_cry.R")
source("./functions/nma_cont.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
source("./functions/nma_utility_functions.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Winbugs stuff
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
model = normal_models()
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
# bugsdir = "C:/Users/TheTimbot/Desktop/WinBUGS14"

# #========================================================================================
# 
# 
# Outcome: Cry Reactivity-----
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
cry_reac_data = NULL

cry_reac_data$pa$gemtc = prep_gem(cry_reac$data)

cry_reac_data$pa$gemtc$data = cry_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(cry_reac$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))



cry_reac_data$pa$network = mtc.network(data.re =cry_reac_data$pa$gemtc$data[,1:4],
                                       studies =cry_reac_data$pa$gemtc$data[,c(1,5:10)])

cry_reac_data$pa$results = mtc.model(cry_reac_data$pa$network, type = "consistency",
                                     linearModel = "random",likelihood = "normal",
                                     link = "identity")
cry_reac_data$pa$results = mtc.run(cry_reac_data$pa$results)

summary(cry_reac_data$pa$results)


cry_reac_data$pa$anohe = mtc.anohe(cry_reac_data$pa$network)
plot(summary(cry_reac_data$pa$anohe))

forest(relative.effect.table(cry_reac_data$pa$results),"drops")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- Informative priors on sigma, huge improvement
# Large diff = 5 points = 1/25= 0.04
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

cry_reac_data$sa1$results = mtc.model(cry_reac_data$pa$network, type = "consistency",
                                      linearModel = "random",likelihood = "normal",
                                      hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.04), 
                                      link = "identity")

cry_reac_data$sa1$results = mtc.run(cry_reac_data$sa1$results)

cry_reac_sucra = sucra(cry_reac_data$sa1$results, direction = -1)
cry_reac_sucra = as.data.frame(cry_reac_sucra) %>% rownames_to_column("treatment")

summary(cry_reac_data$sa1$results)
forest(relative.effect.table(cry_reac_data$sa1$results),"drops")
# save(cry_reac,file = "./cache/cry_reac.rda")

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



source("./analyses/final/heart rate/rop_explore_hr.R")
source("./functions/nma_cont.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
source("./functions/nma_utility_functions.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
model = normal_models()
# bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
bugsdir = "C:/Users/TheTimbot/Desktop/WinBUGS14"



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#no elligible studies for meta analysis
hr_reac_excluded


# hr_reac$bugs$data= prep_wb(data = hr_reac$data,smd = FALSE)
# hr_reac$bugs$pa = nma_cont(hr_reac$bugs$data$wide,hr_reac$bugs$data$wb,hr_reac$bugs$data$treatments,params = params.re, model = list(model$re3,model$re3_inc),
#                             bugsdir = bugsdir, n.iter = 40000, n.burnin = 20000,n.thin = 1, FE = FALSE,debug =F,inc = TRUE)



# save(hr_reac,file = "./cache/hr_reac.rda")

#Informative prior
#Large difference is 5 = 1/25... DIC = 2 so we will use original model


# #========================================================================================
# 
# 
# Outcome: Heart rate Recovery-----
# 
# #========================================================================================


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

hr_recov_data = NULL

hr_recov_data$pa$gemtc = prep_gem(hr_recov$data)

hr_recov_data$pa$gemtc$data = hr_recov_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(hr_recov$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))



hr_recov_data$pa$network = mtc.network(data.re =hr_recov_data$pa$gemtc$data[,1:4],
                                       studies =hr_recov_data$pa$gemtc$data[,c(1,5:10)])



hr_recov_data$pa$results = mtc.model(hr_recov_data$pa$network, type = "consistency",
                                     linearModel = "fixed",likelihood = "normal",
                                     link = "identity")
hr_recov_data$pa$results = mtc.run(hr_recov_data$pa$results)

summary(hr_recov_data$pa$results)


forest(relative.effect.table(hr_recov_data$pa$results),"drops")



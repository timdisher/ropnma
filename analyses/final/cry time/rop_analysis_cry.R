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


source("./analyses/final/cry time/rop_explore_cry.R")

source("./functions/nma_cont.R")
source("./functions/nma_utility_functions.R")
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
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
model = normal_models()


(cry_reac$bugs$data= prep_wb(data = cry_reac$data,smd = FALSE))


cry_reac$bugs$pa = nma_cont(cry_reac$bugs$data$wide,cry_reac$bugs$data$wb,cry_reac$bugs$data$treatments,params = params.re, model = list(model$re2,model$re2_inc),
                           bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 16, FE = FALSE,debug =F,inc = TRUE)






#==== Informative priors
# Large difference = 5 = 1/25 = 0.04, will use for poster but need to revisit

cry_reac$sa1= prep_wb(data = cry_reac$data,smd = FALSE)

cry_reac$sa1$xo = cry_reac$sa1$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>%
  
  left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
  
  mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>%
  
  select(matches("t_"),matches("y_"),matches("se_"),na) %>% select(-y_1)

cry_reac$sa1 = nma_cont(cry_reac$data, cry_reac$sa1$xo, cry_reac$sa1$treatments,params = params.re, model = "re_normal_gaus_2arm_cry_inf.txt",
                            bugsdir = bugsdir, n.iter = 40000, n.burnin = 20000,n.thin = 1, FE = FALSE)



save(cry_reac,file = "./cache/cry_reac.rda")
# #========================================================================================
# 
# 
# Outcome: Cry Recovery-----
# 
# #========================================================================================

#Only Mehta to report narratively
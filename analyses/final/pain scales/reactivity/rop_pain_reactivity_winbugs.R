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

source("./analyses/final/rop_explore_pipp.R")

source("./functions/nma_cont.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_reac$trt_group = fct_infreq(pa_reac$trt_group) #reorders factor from most to least common

pa_reac = droplevels(pa_reac) # Drops factor levels that don't exist (otherwise they are carried over)

#Convert to long format and ensure treatment order is maintained
pa_reac_wb = long_wb(data = pa_reac)


#Correct standard errors for crossovers
pa_reac_wb$wb_xo = pa_reac_wb$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>% left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>% 
  mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>% select(t_1:t_4,y_2:y_4,se_2:se_4,V,na)


#Load models
model = normal_models()


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- Include residual deviance >2
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.fe = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better")
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")

fe.model = nma_cont(pa_reac_wb$wb_xo, pa_reac_wb$treatments,params = params.fe, model = model$fe)

re.model = nma_cont(pa_reac_wb$wb_xo, pa_reac_wb$treatments,params = params.re, model = model$re)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Inconsistency analysis - requires a little extra work in function
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# fe.model.inc = nma_cont_fe(pa_reac_wb$wb_xo, pa_reac_wb$treatments,model = model$fe_inc)
# 
# re.model.inc = nma_cont_fe(pa_reac_wb$wb_xo, pa_reac_wb$treatments,FE = FALSE,model = model$re_inc)
# 
# fe.model.inc
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 ----
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- Include residual deviance >2
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- Include residual deviance >2
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- **Exclude residual deviance >2**
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4 ----
# --------- Include crossovers
# --------- **Standardized mean difference**
# --------- Include imputed means (as SMD)
# --------- Include residual deviance >2
# --------- **Informative priors on sigma**
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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

prep_wb = function(data, smd = FALSE){
data$trt_group = fct_infreq(data$trt_group) %>% droplevels(data)#reorders factor from most to least common

data = droplevels(data) # Drops factor levels that don't exist (otherwise they are carried over)

if(smd == TRUE) long_wb_smd(data) else long_wb(data = data)

}

pa_reac_data = NULL


pa_reac_data$pa = prep_wb(pa_reac)
#Convert to long format and ensure treatment order is maintained



#Correct standard errors for crossovers
pa_reac_data$pa$wb_xo = pa_reac_data$pa$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>% left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>% 
  mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>% select(t_1:t_4,y_2:y_4,se_2:se_4,V,na)


#Load models
model = normal_models()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")

bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"



pa_reac_data$pa$re = nma_cont(pa_reac_data$pa$wb_xo, pa_reac_data$pa$treatments,params = params.re, model = model$re,
                        bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 16, FE = FALSE) 



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 ----
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sens_analysis = function(data,variable,drop,SA = FALSE){
  data = data %>% left_join(rop_data_study[c("studlab","design")],by = "studlab")
  
  if(SA == TRUE){
    data = data %>% filter(data[variable] != drop)} else data = data
    
    pw = pairwise(data = data,treat = trt_group, n = sample_size,mean = mean, sd = std_dev, studlab = studlab, sm = "MD")
    
    nc = netconnection(data = pw, treat1,treat2)
    
    if(nc$n.subnets > 1) return("This sensitivity analysis disconnects the network") else( data = prep_wb(data))
    
    data$xo = data$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>% left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>% 
      mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>% select(t_1:t_4,y_2:y_4,se_2:se_4,V,na)
    nma = nma_cont(data$xo, data$treatments,params = params.re, model = model$re,bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 16, FE = FALSE) 
    
    list(data,nma)
}


pa_reac_data$sa1 = sens_analysis(pa_reac,"design","Crossover") 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa2 = sens_analysis(pa_reac,"imputed_mean","yes")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 ---- Vague priors
# --------- Include crossovers
# --------- Mean difference
# --------- Include imputed means (as SMD)
# --------- include scaled scores
# --------- Vague priors
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa3 = sens_analysis(pa_reac,"scaled_score","yes")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4 ---- Vague priors
# --------- Include crossovers
# --------- **Standardized mean difference**
# --------- Include imputed means (as SMD)
# --------- include scaled scores
# --------- **Informative priors on sigma**
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
pa_reac_data$sa4 = prep_wb(data = pa_reac,smd = TRUE)



pa_reac_data$sa4$re = nma_cont(pa_reac_data$sa4, pa_reac_data$pa$treatments,params = params.re, model = model$re_inf,
                                  bugsdir = bugsdir, n.iter = 100000, n.burnin = 40000,n.thin = 16, FE = FALSE)


data = pa_reac_data$sa4
treatments = pa_reac_data$pa$treatments
params = params.re
model = model$re_inf
bugsdir = bugsdir



# #========================================================================================
# 
# 
# Outcome: Pain Reactivity
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
pa_recov_data = NULL

pa_recov_data$pa$model = sens_analysis(pa_recov)

data = pa_recov
variable = NULL
drop = NULL

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 ----
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_recov_data$sa1 = sens_analysis(pa_recov,"design","Crossover") 

variable = "design"
drop = "Crossover"
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data$sa2 = sens_analysis(pa_recov,"imputed_mean","yes")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 ---- Vague priors
# --------- Include crossovers
# --------- Mean difference
# --------- Include imputed means (as SMD)
# --------- include scaled scores
# --------- Vague priors
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data$sa3 = sens_analysis(pa_recov,"scaled_score","yes")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4 ---- Vague priors
# --------- Include crossovers
# --------- **Standardized mean difference**
# --------- Include imputed means (as SMD)
# --------- include scaled scores
# --------- **Informative priors on sigma**
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
pa_recov_data$sa4 = prep_wb(data = pa_recov,smd = TRUE)



pa_recov_data$sa4$re = nma_cont(pa_recov_data$sa4, pa_recov_data$pa$treatments,params = params.re, model = model$re_inf,
                               bugsdir = bugsdir, n.iter = 100000, n.burnin = 40000,n.thin = 16, FE = FALSE)


data = pa_recov_data$sa4
treatments = pa_recov_data$pa$treatments
params = params.re
model = model$re_inf
bugsdir = bugsdir

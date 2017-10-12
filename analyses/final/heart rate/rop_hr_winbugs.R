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


source("./analyses/final/heart rate/rop_explore_hr.R")

source("./functions/nma_cont.R")
source("./functions/nma_utility_functions.R")
# #========================================================================================
# 
# 
# Outcome: Heart rate Reactivity-----
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
# bugsdir = "C:/Users/TheTimbot/Desktop/WinBUGS14"
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
model = normal_models()

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


hr_recov$bugs$data= prep_wb(data = hr_recov$data,smd = FALSE)


hr_recov$bugs$pa = nma_cont(hr_recov$bugs$data$wide,hr_recov$bugs$data$wb,hr_recov$bugs$data$treatments,params = params.re, model = list(model$re3,model$re3_inc),
                           bugsdir = bugsdir, n.iter = 40000, n.burnin = 20000,n.thin = 1, FE = FALSE,debug =F,inc = TRUE)


#Informative priors - put too much weight on high values

hr_recov$sa1= prep_wb(data = hr_recov$data,smd = FALSE)



hr_recov$sa1 = nma_cont(hr_recov$data, hr_recov$sa1$wb, hr_recov$sa1$treatments,params = params.re, model = "re_normal_gaus_3arm_hr_inf.txt",
                       bugsdir = bugsdir, n.iter = 40000, n.burnin = 20000,n.thin = 1, FE = FALSE)


# save(hr_recov,file = "./cache/hr_recov.rda")

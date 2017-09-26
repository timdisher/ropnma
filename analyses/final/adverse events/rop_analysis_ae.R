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


source("./analyses/final/adverse events/rop_explore_ae.R")

source("./functions/nma_binom.R")
source("./functions/nma_utility_functions.R")
# #========================================================================================
# 
# 
# Outcome: Adverse events Reactivity-----
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
params.fe = c("or", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better")
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
model = normal_models()

ae_reac$bugs$data= prep_wb(data = ae_reac$data_connect,md = FALSE)

ae_reac$bugs$pa = nma_cont(ae_reac$bugs$data$wide,ae_reac$bugs$data$wb,ae_reac$bugs$data$treatments,params = params.fe, model = model$fe,
                           bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 16, FE = FALSE,debug =F,inc = FALSE, cont = FALSE)





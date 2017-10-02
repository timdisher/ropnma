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


source("./analyses/final/spo2/rop_explore_spo2.R")

source("./functions/nma_cont.R")
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
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
model = normal_models()

os_reac$bugs$data= prep_wb(data = os_reac$data,smd = FALSE)


os_reac$bugs$pa = nma_cont(os_reac$bugs$data$wide,os_reac$bugs$data$wb,os_reac$bugs$data$treatments,params = params.re, model = list(model$re2,model$re2_inc),
                           bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 16, FE = FALSE,debug =F,inc = TRUE)






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
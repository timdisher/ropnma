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

prep_wb = function(data, smd = FALSE,convert_contrast = TRUE){
data$trt_group = fct_infreq(data$trt_group) %>% droplevels(data)#reorders factor from most to least common

data = droplevels(data) # Drops factor levels that don't exist (otherwise they are carried over)


if(convert_contrast == TRUE){
if(smd == TRUE) long_wb_smd(data) else long_wb(data = data)}else long_arm_wb(data = data)
  
}

nma = function(data,variable,drop,SA = FALSE, inc = FALSE, models = model){
  data = data %>% left_join(rop_data_study[c("studlab","design")],by = "studlab")
  
  if(SA == TRUE){
    data = data %>% filter(data[variable] != drop)} else data = data
    
    pw = pairwise(data = data,treat = trt_group, n = sample_size,mean = mean, sd = std_dev, studlab = studlab, sm = "MD")
    
    nc = netconnection(data = pw, treat1,treat2)
    
    if(nc$n.subnets > 1) return("This sensitivity analysis disconnects the network") else( data = prep_wb(data))
    
    data$xo = data$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>%
      
      left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
      
      mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>%
      
      select(matches("t_"),matches("y_"),matches("se_"),V,na) %>% select(-y_1)
    
    
    if("y_4" %in% colnames(data$xo))  model = list(models$re,models$re_inc) else model = list(models$re3,models$re3_inc)
    
    if(inc == FALSE){
    nma = nma_cont(data, data$xo, data$treatments,params = params.re, model = model[[1]],bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 10, FE = FALSE, inc = inc)
    
    } else nma = nma_cont(data$wide, data$xo, data$treatments,params = params.re, model = model,bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 10, FE = FALSE, inc = inc)
    list(data = data,nma = nma)
}

treatments = data$treatments
data = data$xo


# #========================================================================================
# 
# 
# Outcome: Pain Reactivity-----
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
model = normal_models()
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"

pa_reac_data = NULL

pa_reac_data$pa = nma(pa_reac, inc = TRUE)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_reac_data$sa1 = nma(pa_reac,"design","Crossover",SA = TRUE) 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa2 = nma(pa_reac,"imputed_mean","yes", SA = TRUE)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 ---- Vague priors
# --------- Include crossovers
# --------- Mean difference
# --------- Include imputed means (as SMD)
# --------- *Exclude scaled scores*
# --------- Vague priors
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa3 = nma(pa_reac,"scaled_score","yes")


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
pa_reac_data$sa4= prep_wb(data = pa_reac,smd = TRUE)



pa_reac_data$sa4 = nma_cont(pa_reac_data$sa4, pa_reac_data$pa$data$treatments,params = params.re, model = model$re_inf,
                                  bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 16, FE = FALSE)

pa_reac_data$sa4





#=========================================
# Sensitivity 5 - Meta-regression on sample size
#=========================================
wf_test = prep_wb(pa_reac)

wf_test$meta = as.data.frame(wf_test$wide %>% rowwise() %>% mutate(x = sum(n_1,n_2,n_3,n_4, na.rm = TRUE)) %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>%
                               
                               left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
                               
                               mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>%
                               
                               select(matches("t_"),matches("y_"),matches("se_"),V,na,x) %>% select(-y_1))


test = nma_winbugs_datalist(wf_test$meta,wf_test$treatments)
test$x = as.vector(wf_test$meta$x)
test$mx = mean(wf_test$meta$x)

params_mr  = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd", "B")

metaregtest = bugs(test,NULL,params_mr,model.file = MODELFILE.re_meta,
                   n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                   bugs.directory = bugsdir, debug = F)

wf_test$meta %>% filter(t_1 == 1) %>% gather(diff,value,y_2:y_4) %>% gather(trt,num, t_2:t_4) %>%
  select(diff,value,num,x) %>% na.omit %>% left_join(test$treatments, by = c("num" = "t")) %>% ggplot(aes(x = x, y = value)) + 
  geom_point() + facet_wrap(~description) + geom_smooth(method = "lm",se = FALSE, colour = "black")

nma_outputs(model = metaregtest,wf_test$treatments)

#=========================================
# Sensitivity 6 - Meta-regression on control arm risk
#=========================================

reac_br_meta = prep_wb(pa_reac)
pa_reac_data$control_risk_data = reac_br_meta$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3),
  se_4 = sd_4/sqrt(n_4)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

(wb_br_meta = nma_winbugs_datalist(test$ready,test$treatments,contrast = FALSE))

wb_br_meta$mx = as.vector(pa_reac_data$control_risk_data %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


metaregtest_armwise = bugs(wb_br_meta,NULL,params_mr,model.file = model$re_arm_meta,
                           n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                           bugs.directory = bugsdir, debug = F)



regress_results = nma_outputs(model = metaregtest_armwise,reac_br_meta$treatments)

test$ready %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                           y3_diff = y_3 - y_1,
                                           y4_diff = y_4 - y_1) %>% gather(diff,value,y2_diff:y4_diff) %>% gather(trt,num, t_2:t_4) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(test$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)






# #========================================================================================
# 
# 
# Outcome: Pain Reactivity-----
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

pa_recov_data$pa = nma(pa_recov, inc = FALSE)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 ----
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_recov_data$sa1 = nma(pa_recov,"design","Crossover", SA = TRUE) 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data$sa2 = nma(pa_recov,"imputed_mean","yes")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 ---- Vague priors
# --------- Include crossovers
# --------- Mean difference
# --------- Include imputed means (as SMD)
# --------- include scaled scores
# --------- Vague priors
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data$sa3 = nma(pa_recov,"scaled_score","yes")


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


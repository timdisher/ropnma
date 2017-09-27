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


source("./analyses/final/pain scales/rop_explore_pipp.R")

source("./functions/nma_cont.R")
source("./functions/nma_utility_functions.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
model = normal_models()
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"



pa_reac_data = NULL


pa_reac_data$pa = nma(pa_reac, inc = TRUE)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1 
# --------- Vague priors
# --------- **Exclude crossovers**
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_reac_data$sa1 = nma(pa_reac,"design","Crossover",SA = TRUE) 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2
# --------- Vague priors
# --------- Include crossovers
# --------- Mean difference outcome
# --------- **Exclude imputed means**
# --------- include scaled scores
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa2 = nma(pa_reac,"imputed_mean","yes", SA = TRUE)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 
# --------- Vague priors
# --------- Include crossovers
# --------- Mean difference
# --------- Include imputed means
# --------- *Exclude scaled scores*
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa3 = nma(pa_reac,"scaled_score","yes", SA = TRUE)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4 ---- 
# --------- Include crossovers
# --------- include scaled scores
# --------- include imputed means
# --------- vague priors
# --------- *exclude posters$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
pa_reac_nopost = pa_reac %>% left_join(rop_data_study[c("studlab","pub_type")], by = "studlab") %>% filter(pub_type == "journal")

pa_reac_data$sa4 = nma(pa_reac_nopost, SA = FALSE)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 5 ---- 
# --------- Include crossovers
# --------- include scaled scores
# --------- include imputed means
# --------- **Informative priors on sigma**
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
pa_reac_data$sa5= prep_wb(data = pa_reac,smd = FALSE)

pa_reac_data$sa5 = nma_cont(pa_reac, pa_reac_data$sa5$wb, pa_reac_data$sa5$treatments,params = params.re, model = "re-normal-gaus_inf_direct.txt",
                                  bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 5, FE = FALSE)


#===============================================================
# Sensitivity 6 - 
# --------- Include crossovers
# --------- include scaled scores
# --------- include imputed means
# --------- vague priors on sigma
# --------- *Weakly informative priors on treatment effects*
# Difference of 3.5 considered to be very large = 1/12.25 = 0.082
#===============================================================
pa_reac_data$sa6= prep_wb(data = pa_reac,smd = FALSE)

pa_reac_data$sa6 = nma_cont(pa_reac, pa_reac_data$sa6$wb, pa_reac_data$sa6$treatments,params = params.re, model = "re-normal-gaus_skep_priors.txt",
                                  bugsdir = bugsdir, n.iter = 200000, n.burnin = 40000,n.thin = 5, FE = FALSE)


#=========================================
# Sensitivity 7 - Meta-regression on high overall risk of bias
#=========================================
pa_reac_data$sa7 = prep_wb(pa_reac)



pa_reac_data$sa7$meta = as.data.frame(pa_reac_data$sa7$wide %>% left_join(rop_data_study %>% select(studlab,oa_rob_sub,design), by = "studlab") %>%
                               
                               left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
                               
                               mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2),
                                      oa_rob_sub = ifelse(oa_rob_sub == "high",1,0)) %>% rename(x = oa_rob_sub) %>%
                               
                               select(matches("t_"),matches("y_"),matches("se_"),V,na,x) %>% select(-y_1))


pa_reac_data$sa7$list = nma_winbugs_datalist(pa_reac_data$sa7$meta,pa_reac_data$sa7$treatments)
pa_reac_data$sa7$list$x = as.vector(pa_reac_data$sa7$meta$x)


params_mr  = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd", "B")

pa_reac_data$sa7$bugs = bugs(pa_reac_data$sa7$list,NULL,params_mr,model.file = model$re_meta,
                   n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                   bugs.directory = bugsdir, debug = F)

pa_reac_data$sa7$meta %>% filter(t_1 == 1) %>% gather(diff,value,y_2:y_4) %>% gather(trt,num, t_2:t_4) %>%
  select(diff,value,num,x) %>% na.omit %>% left_join(wf_test$treatments, by = c("num" = "t")) %>% ggplot(aes(x = x, y = value)) + 
  geom_point() + facet_wrap(~description) + geom_smooth(method = "lm",se = FALSE, colour = "black")

pa_reac_data$sa7$bugs = nma_outputs(model = pa_reac_data$sa7$bugs,pa_reac_data$sa7$treatments)

pa_reac_data$sa7$bugs$B = pa_reac_data$sa7$bugs$bugs[grep("^B$", rownames(pa_reac_data$sa7$bugs$bugs)),]

#=========================================
# Sensitivity 8 - Meta-regression on control arm risk
#=========================================

pa_reac_data$sa8 = prep_wb(pa_reac)

pa_reac_data$sa8$meta_cr = pa_reac_data$sa8$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3),
  se_4 = sd_4/sqrt(n_4)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

(pa_reac_data$sa8$list = nma_winbugs_datalist(pa_reac_data$sa8$meta_cr,pa_reac_data$sa8$treatments,contrast = FALSE))

pa_reac_data$sa8$list$mx = as.vector(pa_reac_data$sa8$meta_cr %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


pa_reac_data$sa8$bugs = bugs(pa_reac_data$sa8$list,NULL,params_mr,model.file = model$re_arm_meta,
                           n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                           bugs.directory = bugsdir, debug = F)



pa_reac_data$sa8$bugs = nma_outputs(model = pa_reac_data$sa8$bugs,pa_reac_data$sa8$treatments)

pa_reac_data$sa8$bugs$B = pa_reac_data$sa8$bugs$bugs[grep("^B$", rownames(pa_reac_data$sa8$bugs$bugs)),]

pa_reac_data$sa8$meta_cr %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                           y3_diff = y_3 - y_1,
                                           y4_diff = y_4 - y_1) %>% gather(diff,value,y2_diff:y4_diff) %>% gather(trt,num, t_2:t_4) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(pa_reac_data$sa8$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)




#=========================================



#==========================================================================
# Power analysis
# Assumptions:
# 1 = I2 is 50%
# 2 = I2 is 70%
# Effect to detect = 2 points (1 MID)
# Assumed SD of PIPP = used Dhaliwhal (largest study)
sd_pipp = sqrt((2.4^2*76+2.1^2*76)/(76+76))

req_samp = round((power.t.test(sig.level = 0.05,power = 0.8,delta = 2, sd = 2.25))$n,0)*2
#==========================================================================

momlinc_netgraph(pa_reac_int,chars$int_char,2)

pa_reac_data$pa$data$wide %>% arrange(n_1)


eff_ss = function(n){
  round(prod(n)/sum(n),0) 
}

power = chars$direct_zeros %>% rename(ctrl = `Treatment Description.x`,
                                      trt = `Treatment Description.y`) %>%  unite(comp,ctrl,trt,sep = " vs ") %>% select(comp, ntot,nstud) 

no_het = power %>% filter(nstud <2) %>% mutate(adj_n = ntot)

het = power %>% filter(nstud >1) %>% mutate(isquared = 0)
  
for(i in 1:length(het$comp)){
  het$isquared[i] = 1-pa_reac_pairwise[names(pa_reac_pairwise) == het$comp[i]][[1]]$I2
}

comb = bind_rows(het,no_het) %>% mutate(isquared = ifelse(is.na(isquared),1,isquared),
                                        adj_n = round(ntot*isquared,0))

comps = comb %>% select(comp,ntot) %>% spread(comp,ntot) #Create wide format of unadjusted ns

comps_adj = comb %>% select(comp,adj_n) %>% spread(comp,adj_n) #Create wide format of adjusted ns

#Adjusted

n_adj = comb[grep("drops vs",comb$comp),c(1,5)] %>%
  mutate(indirect_n_adj =c(
    #drops vs drops sweet
    eff_ss(c(comps_adj$`drops phys vs drops sweet`,comps_adj$`drops vs drops phys`)) +
    eff_ss(c(comps_adj$`drops sweet vs drops sweet mult`,comps_adj$`drops vs drops sweet mult`))+
    eff_ss(c(comps_adj$`drops sweet vs drops.acet`,comps_adj$`drops vs drops.acet`)),

    eff_ss(c(comps_adj$`drops phys vs drops sweet mult`,comps_adj$`drops vs drops phys`)) +
    eff_ss(c(comps_adj$`drops sweet vs drops sweet mult`,comps_adj$`drops vs drops sweet`))+
    eff_ss(c(comps_adj$`drops ebm mult vs drops sweet mult`,comps_adj$`drops vs drops ebm mult`)),
   
   #drops vs drops.acet
   eff_ss(c(comps_adj$`drops sweet vs drops.acet`,comps_adj$`drops vs drops sweet`)),
   
   #drops vs placebo
   0,
   
   #drops vs drops ebm mult
   eff_ss(c(comps_adj$`drops ebm mult vs drops sweet mult`,comps_adj$`drops vs drops sweet mult`)) +
     eff_ss(c(comps_adj$`drops ebm mult vs drops phys`,comps_adj$`drops vs drops phys`)),
   
   #drops vs drops phys
   eff_ss(c(comps_adj$`drops phys vs drops sweet mult`,comps_adj$`drops vs drops sweet mult`)) +
     eff_ss(c(comps_adj$`drops phys vs drops sweet`,comps_adj$`drops vs drops sweet`))+
     eff_ss(c(comps_adj$`drops ebm mult vs drops phys`,comps_adj$`drops vs drops ebm mult`)),
   
   #drops vs drops WFDRI
   0,
   
   #drops vs N20 sweet
   eff_ss(c(comps_adj$`drops sweet vs drops.N2O.sweet`,comps_adj$`drops vs drops sweet`)),
   
   #drops vs sweet
   0,
   
   #drops vs sweet rep
   eff_ss(c(comps_adj$`placebo vs sweet rep`,comps_adj$`drops vs placebo`)),
   
   #drops vs sweet sing
   eff_ss(c(comps_adj$`placebo vs sweet sing`,comps_adj$`drops vs placebo`)))
)

#Unadjusted
n = comb[grep("drops vs",comb$comp),c(1,2)] %>%
  mutate(indirect_n =c(
    #drops vs drops sweet
    eff_ss(c(comps$`drops phys vs drops sweet`,comps$`drops vs drops phys`)) +
      eff_ss(c(comps$`drops sweet vs drops sweet mult`,comps$`drops vs drops sweet mult`))+
      eff_ss(c(comps$`drops sweet vs drops.acet`,comps$`drops vs drops.acet`)),
    
    eff_ss(c(comps$`drops phys vs drops sweet mult`,comps$`drops vs drops phys`)) +
      eff_ss(c(comps$`drops sweet vs drops sweet mult`,comps$`drops vs drops sweet`))+
      eff_ss(c(comps$`drops ebm mult vs drops sweet mult`,comps$`drops vs drops ebm mult`)),
    
    #drops vs drops.acet
    eff_ss(c(comps$`drops sweet vs drops.acet`,comps$`drops vs drops sweet`)),
    
    #drops vs placebo
    0,
    
    #drops vs drops ebm mult
    eff_ss(c(comps$`drops ebm mult vs drops sweet mult`,comps$`drops vs drops sweet mult`)) +
      eff_ss(c(comps$`drops ebm mult vs drops phys`,comps$`drops vs drops phys`)),
    
    #drops vs drops phys
    eff_ss(c(comps$`drops phys vs drops sweet mult`,comps$`drops vs drops sweet mult`)) +
      eff_ss(c(comps$`drops phys vs drops sweet`,comps$`drops vs drops sweet`))+
      eff_ss(c(comps$`drops ebm mult vs drops phys`,comps$`drops vs drops ebm mult`)),
    
    #drops vs drops WFDRI
    0,
    
    #drops vs N20 sweet
    eff_ss(c(comps$`drops sweet vs drops.N2O.sweet`,comps$`drops vs drops sweet`)),
    
    #drops vs sweet
    0,
    
    #drops vs sweet rep
    eff_ss(c(comps$`placebo vs sweet rep`,comps$`drops vs placebo`)),
    
    #drops vs sweet sing
    eff_ss(c(comps$`placebo vs sweet sing`,comps$`drops vs placebo`)))
  )


power_table = left_join(n,n_adj, by = "comp") %>% mutate(tot_eff = ntot + indirect_n,
                                                         tot_eff_adj = adj_n + indirect_n_adj,
                                                         sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                         power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2),
                                                         sample_frac_adj = ifelse(tot_eff_adj/req_samp >1,">100",round(tot_eff_adj/req_samp*100,2)),
                                                         power_adj = round((power.t.test(tot_eff_adj/2,delta = 2, sd = 2.25))$power,2)) %>% 
  
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,adj_n,indirect_n_adj,tot_eff_adj,sample_frac_adj,power_adj)



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

pa_recov_data$pa = nma(pa_recov, inc = TRUE)


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


#=========================================
# Sensitivity 5 - Meta-regression on sample size
#=========================================
pa_recov_data$sa5$data = prep_wb(pa_recov)

pa_recov_data$sa5$data$meta = as.data.frame(pa_recov_data$sa5$data$wide %>% rowwise() %>% mutate(x = sum(n_1,n_2,n_3, na.rm = TRUE)) %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>%
                               
                               left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
                               
                               mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>%
                               
                               select(matches("t_"),matches("y_"),matches("se_"),V,na,x) %>% select(-y_1))


pa_recovsa5_data = nma_winbugs_datalist(pa_recov_data$sa5$data$meta,pa_recov_data$sa5$data$treatments)
pa_recovsa5_data$x = as.vector(pa_recov_data$sa5$data$meta$x)
pa_recovsa5_data$mx = mean(pa_recovsa5_data$x)

params_mr  = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd", "B")

recov_metareg_samplesize = bugs(pa_recovsa5_data,NULL,params_mr,model.file = model$re3_meta,
                   n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                   bugs.directory = bugsdir, debug = F)


pa_recov_data$sa5$tables = nma_outputs(model = recov_metareg_samplesize,pa_recov_data$sa5$data$treatments)
pa_recov_data$sa5$tables$B = as.data.frame(recov_metareg_samplesize$summary) %>% rownames_to_column("parameter") %>% filter(parameter == "B")

pa_recov_data$sa5$data$meta %>% filter(t_1 == 1) %>% select(-t_1) %>% gather(diff,value,matches("y_")) %>% gather(trt,num, matches("t_")) %>%
  select(diff,value,num,x) %>% na.omit %>% left_join(pa_recov_data$sa5$data$treatments, by = c("num" = "t")) %>% ggplot(aes(x = x, y = value)) + 
  geom_point() + facet_wrap(~description) + geom_smooth(method = "lm",se = FALSE, colour = "black")

#=========================================
# Sensitivity 6 - Meta-regression on control arm risk
#=========================================

pa_recov_data$sa6$data = prep_wb(pa_recov)

pa_recov_data$sa6$data$meta = pa_recov_data$sa6$data$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

pa_recovsa6_data = nma_winbugs_datalist(pa_recov_data$sa6$data$meta,
                                        pa_recov_data$sa6$data$treatments,contrast = FALSE)

pa_recovsa6_data$mx = as.vector(pa_recov_data$sa6$data$meta %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


recov_metareg_baselinerisk = bugs(pa_recovsa6_data,NULL,params_mr,model.file = model$re_arm_meta,
                           n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                           bugs.directory = bugsdir, debug = F)



pa_recov_data$sa6$tables = nma_outputs(model = recov_metareg_baselinerisk,pa_recov_data$sa6$data$treatments)

pa_recov_data$sa6$data$meta %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                           y3_diff = y_3 - y_1) %>% gather(diff,value,ends_with("_diff")) %>% select(-t_1) %>% gather(trt,num, matches("t_")) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(pa_recov_data$sa6$data$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)


#==========================================================================
# Power analysis
# Assumptions:
# 1 = I2 is 50%
# 2 = I2 is 70%
# Effect to detect = 2 points (1 MID)
# Assumed SD of PIPP = used Dhaliwhal (largest study)

#==========================================================================
momlinc_netgraph(pa_recov_int,recov_chars$int_char,2)



recov_power = recov_chars$direct_zeros %>% rename(ctrl = `Treatment Description.x`,
                                      trt = `Treatment Description.y`) %>%  unite(comp,ctrl,trt,sep = " vs ") %>% select(comp, ntot,nstud) 

recov_no_het = recov_power %>% filter(nstud <2) %>% mutate(adj_n = ntot)

recov_het = recov_power %>% filter(nstud >1) %>% mutate(isquared = 0)

for(i in 1:length(recov_het$comp)){
  recov_het$isquared[i] = 1-pa_recov_pairwise[names(pa_recov_pairwise) == recov_het$comp[i]][[1]]$I2
}

recov_comb = bind_rows(recov_het,recov_no_het) %>% mutate(isquared = ifelse(is.na(isquared),1,isquared),
                                        adj_n = round(ntot*isquared,0))

recov_comps = recov_comb %>% select(comp,ntot) %>% spread(comp,ntot) #Create wide format of unadjusted ns

recov_comps_adj = recov_comb %>% select(comp,adj_n) %>% spread(comp,adj_n) #Create wide format of adjusted ns

#Adjusted
recov_n = recov_comb[grep("drops vs",recov_comb$comp),c(1,5)] %>%
  mutate(indirect_n_adj =c(
    #drops vs drops sweet
    eff_ss(c(recov_comps_adj$`drops sweet vs drops.acet`,recov_comps_adj$`drops vs drops.acet`)),
    
    #drops vs drops acet
    eff_ss(c(recov_comps_adj$`drops sweet vs drops.acet`,recov_comps_adj$`drops vs drops sweet`)),
    
    #drops vs drops ebm mult
    eff_ss(c(recov_comps_adj$`drops ebm mult vs drops sweet mult`,recov_comps_adj$`drops vs drops sweet mult`)),
    
    #drops vs drops morph
    eff_ss(c(recov_comps_adj$`drops morph vs drops.acet`,recov_comps_adj$`drops vs drops.acet`)),
    
    #drops vs drops phys
    eff_ss(c(recov_comps_adj$`drops phys vs drops sweet mult`,recov_comps_adj$`drops vs drops sweet mult`)) +
      eff_ss(c(recov_comps_adj$`drops ebm mult vs drops phys`,recov_comps_adj$`drops vs drops ebm mult`)),
    
    #drops vs drops sweet mult
    eff_ss(c(recov_comps_adj$`drops ebm mult vs drops sweet mult`,recov_comps_adj$`drops vs drops ebm mult`)),
    
    #drops vs phys
    0,
    
    #drops vs placebo
    0,
    
    #drops vs sweet
    0
    
  ))


#Unadjusted

recov_n = recov_comb[grep("drops vs",recov_comb$comp),c(1,5)] %>%
  mutate(indirect_n =c(
    #drops vs drops sweet
    eff_ss(c(recov_comps$`drops sweet vs drops.acet`,recov_comps$`drops vs drops.acet`)),
    
    #drops vs drops acet
    eff_ss(c(recov_comps$`drops sweet vs drops.acet`,recov_comps$`drops vs drops sweet`)),
    
    #drops vs drops ebm mult
    eff_ss(c(recov_comps$`drops ebm mult vs drops sweet mult`,recov_comps$`drops vs drops sweet mult`)),
    
    #drops vs drops morph
    eff_ss(c(recov_comps$`drops morph vs drops.acet`,recov_comps$`drops vs drops.acet`)),
    
    
    #drops vs drops phys
    eff_ss(c(recov_comps$`drops phys vs drops sweet mult`,recov_comps$`drops vs drops sweet mult`)) +
      eff_ss(c(recov_comps$`drops ebm mult vs drops phys`,recov_comps$`drops vs drops ebm mult`)),
    
    #drops vs drops sweet mult
    eff_ss(c(recov_comps$`drops ebm mult vs drops sweet mult`,recov_comps$`drops vs drops ebm mult`)),
    
    #drops vs phys
    0,
    
    #drops vs placebo
    0,
    
    #drops vs sweet
    0
    
  ))

recov_power_table = left_join(recov_n,recov_n_adj, by = "comp") %>% mutate(tot_eff = ntot + indirect_n,
                                                         tot_eff_adj = adj_n + indirect_n_adj,
                                                         sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                         power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2),
                                                         sample_frac_adj = ifelse(tot_eff_adj/req_samp >1,">100",round(tot_eff_adj/req_samp*100,2)),
                                                         power_adj = round((power.t.test(tot_eff_adj/2,delta = 2, sd = 2.25))$power,2)) %>% 
  
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,adj_n,indirect_n_adj,tot_eff_adj,sample_frac_adj,power_adj)

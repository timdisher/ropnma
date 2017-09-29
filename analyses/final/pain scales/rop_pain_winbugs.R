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
# bugsdir = "C:/Users/TheTimbot/Desktop/WinBUGS14"


pa_reac_data = NULL


pa_reac_data$pa = nma(pa_reac, inc = FALSE)

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
pa_reac_data$pa$data$wide %>% arrange(n_1)
sd_pipp = sqrt((2.4^2*75+2.1^2*75)/(75+75))
req_samp = round((power.t.test(sig.level = 0.05,power = 0.8,delta = 2, sd = 2.25))$n,0)*2
#==========================================================================

momlinc_netgraph(pa_reac_int,chars$int_char,2)



eff_ss = function(n){
  round(prod(n)/sum(n),0) 
}

power = chars$direct_zeros %>% rename(ctrl = `Treatment Description.x`,
                                      trt = `Treatment Description.y`) %>%  unite(comp,trt,ctrl,sep = " vs ") %>% select(comp, ntot,nstud) %>% mutate(num = seq(1,66,1))

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

n_adj = comb[grep("vs drops$",comb$comp),c(1,4,6)] %>%
  mutate(indirect_n_adj =c(
    #drops sweet vs drops
    eff_ss(c(comps_adj$`drops sweet vs drops phys`,comps_adj$`drops phys vs drops`)) +
    eff_ss(c(comps_adj$`drops sweet vs drops sweet mult`,comps_adj$`drops sweet mult vs drops`))+
    eff_ss(c(comps_adj$`drops.acet vs drops sweet`,comps_adj$`drops.acet vs drops`)),

    eff_ss(c(comps_adj$`drops sweet vs drops phys`,comps_adj$`drops phys vs drops`)) +
    eff_ss(c(comps_adj$`drops sweet vs drops sweet mult`,comps_adj$`drops sweet vs drops`))+
    eff_ss(c(comps_adj$`drops ebm mult vs drops sweet mult`,comps_adj$`drops ebm mult vs drops`)),
   
   #drops.acet vs drops
   eff_ss(c(comps_adj$`drops.acet vs drops sweet`,comps_adj$`drops sweet vs drops`)),
   
   #placebo vs drops
   0,
   
   #drops ebm mult vs drops
   eff_ss(c(comps_adj$`drops ebm mult vs drops sweet mult`,comps_adj$`drops sweet mult vs drops`)) +
     eff_ss(c(comps_adj$`drops ebm mult vs drops phys`,comps_adj$`drops phys vs drops`)),
   
   #drops phys vs drops
   eff_ss(c(comps_adj$`drops sweet vs drops phys`,comps_adj$`drops sweet mult vs drops`)) +
     eff_ss(c(comps_adj$`drops sweet vs drops phys`,comps_adj$`drops sweet vs drops`))+
     eff_ss(c(comps_adj$`drops ebm mult vs drops phys`,comps_adj$`drops ebm mult vs drops`)),
   
   #drops vs drops WFDRI
   0,
   
   #drops vs N20 sweet
   eff_ss(c(comps_adj$`drops.N2O.sweet vs drops sweet`,comps_adj$`drops sweet vs drops`)),
   
   #drops vs sweet
   0,
   
   #drops vs sweet rep
   eff_ss(c(comps_adj$`sweet rep vs placebo`,comps_adj$`placebo vs drops`)),
   
   #drops vs sweet sing
   eff_ss(c(comps_adj$`sweet sing vs placebo`,comps_adj$`placebo vs drops`)))
)

#Unadjusted
n = comb[grep("vs drops$",comb$comp),c(1,2,4)] %>%
  mutate(indirect_n =c(
    #drops sweet vs drops
    eff_ss(c(comps$`drops sweet vs drops phys`,comps$`drops phys vs drops`)) +
      eff_ss(c(comps$`drops sweet vs drops sweet mult`,comps$`drops sweet mult vs drops`))+
      eff_ss(c(comps$`drops.acet vs drops sweet`,comps$`drops.acet vs drops`)),
    
    eff_ss(c(comps$`drops sweet vs drops phys`,comps$`drops phys vs drops`)) +
      eff_ss(c(comps$`drops sweet vs drops sweet mult`,comps$`drops sweet vs drops`))+
      eff_ss(c(comps$`drops ebm mult vs drops sweet mult`,comps$`drops ebm mult vs drops`)),
    
    #drops.acet vs drops
    eff_ss(c(comps$`drops.acet vs drops sweet`,comps$`drops sweet vs drops`)),
    
    #placebo vs drops
    0,
    
    #drops ebm mult vs drops
    eff_ss(c(comps$`drops ebm mult vs drops sweet mult`,comps$`drops sweet mult vs drops`)) +
      eff_ss(c(comps$`drops ebm mult vs drops phys`,comps$`drops phys vs drops`)),
    
    #drops phys vs drops
    eff_ss(c(comps$`drops sweet vs drops phys`,comps$`drops sweet mult vs drops`)) +
      eff_ss(c(comps$`drops sweet vs drops phys`,comps$`drops sweet vs drops`))+
      eff_ss(c(comps$`drops ebm mult vs drops phys`,comps$`drops ebm mult vs drops`)),
    
    #drops vs drops WFDRI
    0,
    
    #drops vs N20 sweet
    eff_ss(c(comps$`drops.N2O.sweet vs drops sweet`,comps$`drops sweet vs drops`)),
    
    #drops vs sweet
    0,
    
    #drops vs sweet rep
    eff_ss(c(comps$`sweet rep vs placebo`,comps$`placebo vs drops`)),
    
    #drops vs sweet sing
    eff_ss(c(comps$`sweet sing vs placebo`,comps$`placebo vs drops`)))
  )





power_table = left_join(n,n_adj, by = c("comp","num")) %>% mutate(tot_eff = ntot + indirect_n,
                                                         tot_eff_adj = adj_n + indirect_n_adj,
                                                         sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                         power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2),
                                                         sample_frac_adj = ifelse(tot_eff_adj/req_samp >1,">100",round(tot_eff_adj/req_samp*100,2)),
                                                         power_adj = round((power.t.test(tot_eff_adj/2,delta = 2, sd = 2.25))$power,2)) %>% 
  
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,adj_n,indirect_n_adj,tot_eff_adj,sample_frac_adj,power_adj,num)

retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000) { 
  z <- qt(1-alpha/2, df)
  p.hi <- 1 - pt(z-A/s, df)
  p.lo <- pt(-z-A/s, df)  
  power <- p.hi + p.lo
  typeS <- p.lo/power
  estimate <- A + s*rt(n.sims,df)
  significant <- abs(estimate) > s*z
  exaggeration <- mean(abs(estimate)[significant])/A
  return(list(power=power,
              typeS=typeS, exaggeration=exaggeration))
}


power_table = power_table %>% mutate(post_sd = pa_reac_data$pa$nma$bugs[1:11,2],
                                    gelman_power = 0,
                                    gelman_n = 0,
                                    gelman_m = 0)

for(i in seq_along(power_table$comp)){
  power_table$gelman_power[i] = retrodesign(2,power_table$post_sd[i])[["power"]]
  power_table$gelman_m[i] = retrodesign(2,power_table$post_sd[i])[["exaggeration"]]
  power_table$gelman_n[i] = (power.t.test(power = power_table$gelman_power[i],delta = 2, sd = 2.25))$n*2
  
}



# Forest plot vs ref======
power_table = power_table %>% arrange(num) 
pa_reac_forest_data = as.data.frame(pa_reac_data$pa$nma$comp[1:11,1:3])


pa_reac_forest_data = bind_cols(power_table,pa_reac_forest_data) %>% select(-one_of(c("sample_frac","tot_eff","sample_frac_adj","adj_n",
                                                                                      "indirect_n_adj","power","num",
                                                                                      "Comparison (Trt A vs. Trt B)","tot_eff_adj",
                                                                                      "power_adj"))) %>%
  mutate(comp = c("Drops + sweet taste mult","Drops + phys","Drops + sweet taste", "Placebo",
                  "Drops + ebm mult","Drops + acetaminophen", "Drops + WFDRI",
                  "Drops + N2O + sweet taste","Sweet taste alone","Repeated sweet taste",
                  "Sweet tate + singing")) %>% arrange(as.numeric(as.character(`Mean Difference of Trt A vs. Trt B`))) 

pa_reac_plot_data = pa_reac_forest_data[,c(8,9)] %>%separate(`95% CrI of Mean Difference`,c("lower","upper"),sep = " to ") %>% rename(mean = `Mean Difference of Trt A vs. Trt B`)
pa_reac_table_data = pa_reac_forest_data %>%  mutate(cri = paste("(",`95% CrI of Mean Difference`,")",sep="")) %>% select(-`95% CrI of Mean Difference`) %>% unite(mean_cri,`Mean Difference of Trt A vs. Trt B`,cri,sep = " ") %>%
  select(-post_sd) 

pa_reac_table_data$gelman_n = round(pa_reac_table_data$gelman_n,0)
pa_reac_table_data$gelman_m = round(pa_reac_table_data$gelman_m,2)
pa_reac_table_data$gelman_power = round(pa_reac_table_data$gelman_power,2)

pa_reac_table_data = pa_reac_table_data %>%  mutate(sig = c("yes","no","yes",rep("no",8)),
         gelman_m = ifelse(sig == "yes",gelman_m,"NA"))  %>% select(comp,ntot,indirect_n,gelman_n,gelman_power,gelman_m,mean_cri)



pdf("./figs/final/pain scales reactivity/nma_pa_reac_power_forest.pdf", onefile = FALSE, width = 12, height = 5)
forestplot(rbind(c("Comparison","Direct","Effective 
indirect
","Heterogeneity
adjusted N
                   ","Power","Exaggeration 
factor
                   ","Mean Difference (95% CrI)"),pa_reac_table_data),
           c(NA,as.numeric(as.character(pa_reac_plot_data$mean))),
           c(NA,as.numeric(as.character(pa_reac_plot_data$lower))),
           c(NA,as.numeric(as.character(pa_reac_plot_data$upper))),
           graph.pos = 7, graphwidth = unit(50,'mm'),
           is.summary = c(TRUE,rep(FALSE,11)),
           hrzl_lines = gpar(col="#444444"),
           align = c("l",rep("c",5),"l"),
           boxsize = 0.4,
           colgap = unit(3,"mm"),
           lineheight = unit(7,"mm"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 11, family = "calibri")),
           col = fpColors(box = "mediumpurple", line = "midnightblue"),
           xlab = "PIPP Score",
           title = "PIPP Reactivity (Intervention vs Anesthetic eye drops alone)"
)
dev.off()




# #========================================================================================
# 
# 
# Outcome: Pain Recovery-----
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
                                      trt = `Treatment Description.y`) %>%  unite(comp,trt,ctrl,sep = " vs ") %>% select(comp, ntot,nstud) %>% mutate(num = seq(1,45,1))

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
recov_n_adj = recov_comb[grep("vs drops$",recov_comb$comp),c(1,4,6)] %>%
  mutate(indirect_n_adj =c(
    #drops sweet vs drops
    eff_ss(c(recov_comps_adj$`drops.acet vs drops sweet`,recov_comps_adj$`drops.acet vs drops`)),
    
    #drops vs drops acet
    eff_ss(c(recov_comps_adj$`drops.acet vs drops sweet`,recov_comps_adj$`drops sweet vs drops`)),
    
    #drops ebm mult vs drops
    eff_ss(c(recov_comps_adj$`drops ebm mult vs drops sweet mult`,recov_comps_adj$`drops sweet mult vs drops`)),
    
    #drops vs drops morph
    eff_ss(c(recov_comps_adj$`drops morph vs drops.acet`,recov_comps_adj$`drops.acet vs drops`)),
    
    #drops phys vs drops
    eff_ss(c(recov_comps_adj$`drops sweet vs drops phys`,recov_comps_adj$`drops sweet mult vs drops`)) +
      eff_ss(c(recov_comps_adj$`drops ebm mult vs drops phys`,recov_comps_adj$`drops ebm mult vs drops`)),
    
    #drops sweet mult vs drops
    eff_ss(c(recov_comps_adj$`drops ebm mult vs drops sweet mult`,recov_comps_adj$`drops ebm mult vs drops`)),
    
    #drops vs phys
    0,
    
    #placebo vs drops
    0,
    
    #drops vs sweet
    0
    
  ))


#Unadjusted

recov_n = recov_comb[grep("vs drops$",recov_comb$comp),c(1,2,4)] %>%
  mutate(indirect_n =c(
    #drops sweet vs drops
    eff_ss(c(recov_comps$`drops.acet vs drops sweet`,recov_comps$`drops.acet vs drops`)),
    
    #drops vs drops acet
    eff_ss(c(recov_comps$`drops.acet vs drops sweet`,recov_comps$`drops sweet vs drops`)),
    
    #drops ebm mult vs drops
    eff_ss(c(recov_comps$`drops ebm mult vs drops sweet mult`,recov_comps$`drops sweet mult vs drops`)),
    
    #drops vs drops morph
    eff_ss(c(recov_comps$`drops morph vs drops.acet`,recov_comps$`drops.acet vs drops`)),
    
    
    #drops phys vs drops
    eff_ss(c(recov_comps$`drops sweet vs drops phys`,recov_comps$`drops sweet mult vs drops`)) +
      eff_ss(c(recov_comps$`drops ebm mult vs drops phys`,recov_comps$`drops ebm mult vs drops`)),
    
    #drops sweet mult vs drops
    eff_ss(c(recov_comps$`drops ebm mult vs drops sweet mult`,recov_comps$`drops ebm mult vs drops`)),
    
    #drops vs phys
    0,
    
    #placebo vs drops
    0,
    
    #drops vs sweet
    0
    
  ))

recov_power_table = left_join(recov_n,recov_n_adj, by = c("comp","num")) %>% mutate(tot_eff = ntot + indirect_n,
                                                         tot_eff_adj = adj_n + indirect_n_adj,
                                                         sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                         power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2),
                                                         sample_frac_adj = ifelse(tot_eff_adj/req_samp >1,">100",round(tot_eff_adj/req_samp*100,2)),
                                                         power_adj = round((power.t.test(tot_eff_adj/2,delta = 2, sd = 2.25))$power,2)) %>% 
  
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,adj_n,indirect_n_adj,tot_eff_adj,sample_frac_adj,power_adj,num)


recov_power_table = recov_power_table %>% mutate(post_sd = pa_recov_data$pa$nma$bugs[1:9,2],
                                     gelman_power = 0,
                                     gelman_n = 0,
                                     gelman_m = 0)

for(i in seq_along(recov_power_table$comp)){
  recov_power_table$gelman_power[i] = retrodesign(2,recov_power_table$post_sd[i])[["power"]]
  recov_power_table$gelman_m[i] = retrodesign(2,recov_power_table$post_sd[i])[["exaggeration"]]
  recov_power_table$gelman_n[i] = (power.t.test(power = recov_power_table$gelman_power[i],delta = 2, sd = 2.25))$n*2
  
}

# Forest plot vs ref======
recov_power_table = recov_power_table %>% arrange(num) 
pa_recov_forest_data = as.data.frame(pa_recov_data$pa$nma$comp[1:9,1:3])


pa_recov_forest_data = bind_cols(recov_power_table,pa_recov_forest_data) %>% select(-one_of(c("sample_frac","tot_eff","sample_frac_adj","adj_n",
                                                                                      "indirect_n_adj","power","num",
                                                                                      "Comparison (Trt A vs. Trt B)","tot_eff_adj",
                                                                                      "power_adj"))) %>%
  mutate(comp = c("Drops + phys","Drops + sweet taste mult","Drops + sweet taste","Drops + ebm mult",
                  "Drops + acetaminophen","Drops + morphine","NNS alone","Placebo", 
                  "Sweet taste alone")) %>% arrange(as.numeric(as.character(`Mean Difference of Trt A vs. Trt B`))) 

pa_recov_plot_data = pa_recov_forest_data[,c(8,9)] %>%separate(`95% CrI of Mean Difference`,c("lower","upper"),sep = " to ") %>% rename(mean = `Mean Difference of Trt A vs. Trt B`)
pa_recov_table_data = pa_recov_forest_data %>%  mutate(cri = paste("(",`95% CrI of Mean Difference`,")",sep="")) %>% select(-`95% CrI of Mean Difference`) %>% unite(mean_cri,`Mean Difference of Trt A vs. Trt B`,cri,sep = " ") %>%
  select(-post_sd) 

pa_recov_table_data$gelman_n = round(pa_recov_table_data$gelman_n,0)
pa_recov_table_data$gelman_m = round(pa_recov_table_data$gelman_m,2)
pa_recov_table_data$gelman_power = round(pa_recov_table_data$gelman_power,2)

pa_recov_table_data = pa_recov_table_data %>%  mutate(sig = c("yes","yes",rep("no",7)),
                                                    gelman_m = ifelse(sig == "yes",gelman_m,"NA"))  %>% select(comp,ntot,indirect_n,gelman_n,gelman_power,gelman_m,mean_cri)



pdf("./figs/final/pain scales recovery/nma_pa_recov_power_forest.pdf", onefile = FALSE, width = 12, height = 5)
forestplot(rbind(c("Comparison","Direct","Effective 
indirect
","Heterogeneity
adjusted N
                   ","Power","Exaggeration 
factor
                   ","Mean Difference (95% CrI)"),pa_recov_table_data),
           c(NA,as.numeric(as.character(pa_recov_plot_data$mean))),
           c(NA,as.numeric(as.character(pa_recov_plot_data$lower))),
           c(NA,as.numeric(as.character(pa_recov_plot_data$upper))),
           graph.pos = 7, graphwidth = unit(50,'mm'),
           is.summary = c(TRUE,rep(FALSE,11)),
           hrzl_lines = gpar(col="#444444"),
           align = c("l",rep("c",5),"l"),
           boxsize = 0.4,
           colgap = unit(3,"mm"),
           lineheight = unit(7,"mm"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 11, family = "calibri")),
           col = fpColors(box = "mediumpurple", line = "midnightblue"),
           xlab = "PIPP Score",
           title = "PIPP Recovery (Intervention vs Anesthetic eye drops alone)"
)
dev.off()

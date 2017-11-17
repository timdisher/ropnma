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
library(tidyverse)
library(netmeta)
library(stargazer)
library(reshape2)
library(forcats)
library(scales)
library(forestplot)
library(gemtc)
library(reshape2)
library(R2jags)



source("./analyses/final/pain scales/rop_explore_pipp.R")
source("./functions/nma_cont.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
source("./functions/nma_utility_functions.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
params.re = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd")
model = normal_models()
bugsdir = "C:/Users/dishtc/Desktop/WinBUGS14"
# bugsdir = "C:/Users/TheTimbot/Desktop/WinBUGS14"


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

pa_reac_data = NULL

pa_reac_data$pa$gemtc = prep_gem(pa_reac)


pa_reac_data$pa$gemtc$data = pa_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(pa_reac[c("studlab","imputed_mean","scaled_score")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "high",1,0),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1))



pa_reac_data$pa$network = mtc.network(data.re =pa_reac_data$pa$gemtc$data[,1:4],
                                      studies =pa_reac_data$pa$gemtc$data[,c(1,5:9)])

pa_reac_data$pa$results = mtc.model(pa_reac_data$pa$network, type = "consistency",
                                    linearModel = "random",likelihood = "normal",
                                    link = "identity")
pa_reac_data$pa$results = mtc.run(pa_reac_data$pa$results)

summary(pa_reac_data$pa$results)

pa_reac_data$pa$anohe = mtc.anohe(pa_reac_data$pa$network)
# plot(summary(pa_reac_data$pa$anohe))
forest(relative.effect.table(pa_reac_data$pa$results),"drops")

pa_reac_sucra = sucra(pa_reac_data$pa$results, direction = -1)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- informative priors on sigma, used 4 point difference as very large
# ----- Large improvement in sd, small difference in DIC. Use this moving forward.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa1$results = mtc.model(pa_reac_data$pa$network, type = "consistency",
                                     linearModel = "random",likelihood = "normal",
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625), 
                                     link = "identity")

pa_reac_data$sa1$results = mtc.run(pa_reac_data$sa1$results)

summary(pa_reac_data$sa1$results)


forest(relative.effect.table(pa_reac_data$sa1$results),"drops") 

pa_reac_sucra = sucra(pa_reac_data$sa1$results, direction = -1)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 
# ----- Meta-regression on cross-over design
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sa2_regressor = list(coefficient = "shared",
                     variable = "design",
                     control = "drops")

pa_reac_data$sa2$results = mtc.model(pa_reac_data$pa$network, type = "regression",
                                    linearModel = "random",likelihood = "normal",
                                    regressor = sa2_regressor,
                                    hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625),
                                    link = "identity")


pa_reac_data$sa2$results = mtc.run(pa_reac_data$sa2$results)

summary(pa_reac_data$sa2$results)

forest(relative.effect.table(pa_reac_data$sa2$results),"drops")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 
# ----- Meta-regression on imputed mean
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


sa3_regressor = list(coefficient = "shared",
                     variable = "imputed_mean",
                     control = "drops")



pa_reac_data$sa3$results = mtc.model(pa_reac_data$pa$network, type = "regression",
                                     linearModel = "random",likelihood = "normal",
                                     regressor = sa3_regressor,
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625),
                                     link = "identity")


pa_reac_data$sa3$results = mtc.run(pa_reac_data$sa3$results)

summary(pa_reac_data$sa3$results)

forest(relative.effect.table(pa_reac_data$sa3$results),"drops")

# plot(pa_reac_data$sa3$results)

#
# plot(mtc.deviance(pa_reac_data$sa3$results)$dev.re)
# 
# plot(mtc.deviance(pa_reac_data$pa$results)$dev.re)
#=------Drop imputed mean----
# pa_reac_data$sa3_drop$gemtc = prep_gem(pa_reac %>% filter(imputed_mean == "no"))
# 
# pa_reac_data$sa3_drop$network = mtc.network(data.re =pa_reac_data$sa3_drop$gemtc$data)
# 
# 
# pa_reac_data$sa3_drop$results = mtc.model(pa_reac_data$sa3_drop$network, type = "consistency",
#                                      linearModel = "random",likelihood = "normal",
#                                      hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.625),
#                                      link = "identity")
# 
# pa_reac_data$sa3_drop$results = mtc.run(pa_reac_data$sa3_drop$results)
# 
# 
# summary(pa_reac_data$sa3_drop$results)
# forest(relative.effect.table(pa_reac_data$sa3_drop$results),"drops")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4
# ----- Meta-regression on scaled scores (not enough data, done removing instead)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data$sa4$gemtc = prep_gem(pa_reac %>% filter(scaled_score == "no"))

pa_reac_data$sa4$network = mtc.network(data.re =pa_reac_data$sa4$gemtc$data)
                                       

pa_reac_data$sa4$results = mtc.model(pa_reac_data$sa4$network, type = "consistency",
                                     linearModel = "random",likelihood = "normal",
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625),
                                     link = "identity")

pa_reac_data$sa4$results = mtc.run(pa_reac_data$sa4$results)

                                      
summary(pa_reac_data$sa4$results)
forest(relative.effect.table(pa_reac_data$sa4$results),"drops")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 5
# ----- Meta-regression on pub type (not enough data, done removing instead)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
pa_reac_data$sa5$gemtc = prep_gem(pa_reac %>% left_join(rop_data_study[c("studlab","pub_type")],by = "studlab") %>% filter(pub_type == "journal"))

pa_reac_data$sa5$network = mtc.network(data.re =pa_reac_data$sa5$gemtc$data)


pa_reac_data$sa5$results = mtc.model(pa_reac_data$sa5$network, type = "consistency",
                                     linearModel = "random",likelihood = "normal",
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625),
                                     link = "identity")

pa_reac_data$sa5$results = mtc.run(pa_reac_data$sa5$results)


summary(pa_reac_data$sa5$results)


#=========================================
# Sensitivity 6 - Meta-regression on high overall risk of bias
#=========================================

sa6_regressor = list(coefficient = "shared",
                     variable = "oa_rob_sub",
                     control = "drops")



pa_reac_data$sa6$results = mtc.model(pa_reac_data$pa$network, type = "regression",
                                     linearModel = "random",likelihood = "normal",
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.625),
                                     regressor = sa6_regressor,
                                     link = "identity")

pa_reac_data$sa6$results = mtc.run(pa_reac_data$sa6$results)

summary(pa_reac_data$sa6$results)


#=========================================
# Sensitivity 7 - Meta-regression on control arm risk
#=========================================
params_mr  = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better","sd", "B")


pa_reac_data$sa7 = prep_wb(pa_reac)

pa_reac_data$sa7$meta_cr = pa_reac_data$sa7$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3),
  se_4 = sd_4/sqrt(n_4)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

(pa_reac_data$sa7$list = nma_winbugs_datalist(pa_reac_data$sa7$meta_cr,pa_reac_data$sa7$treatments,contrast = FALSE))

pa_reac_data$sa7$list$mx = as.vector(pa_reac_data$sa7$meta_cr %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


pa_reac_data$sa7$bugs = jags(pa_reac_data$sa7$list, NULL, params_mr,model.file = "./jags_models/re_normal_armdata_meta_jags.txt",
                             n.chains = 3, n.iter = 40000, n.burnin = 20000, n.thin = 1)



pa_reac_data$sa7$bugs = nma_outputs(model = pa_reac_data$sa7$bugs$BUGSoutput,pa_reac_data$sa7$treatments)

pa_reac_data$sa7$bugs$B = pa_reac_data$sa7$bugs$bugs[grep("^B$", rownames(pa_reac_data$sa7$bugs$bugs)),]

pa_reac_data$sa7$meta_cr %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                                         y3_diff = y_3 - y_1,
                                                         y4_diff = y_4 - y_1) %>% gather(diff,value,y2_diff:y4_diff) %>% gather(trt,num, t_2:t_4) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(pa_reac_data$sa7$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)



#==========================================================================
# Power analysis
# Assumptions:
# 1 = I2 is 50%
# 2 = I2 is 70%
# Effect to detect = 2 points (1 MID)
# Assumed SD of PIPP = used Dhaliwhal (largest study)
pa_reac_data$pa$gemtc$input %>% arrange(-n)
sd_pipp = sqrt((2.4^2*75+2.1^2*75)/(75+75))
req_samp = round((power.t.test(sig.level = 0.05,power = 0.8,delta = 2, sd = 2.25))$n,0)*2
#==========================================================================


eff_ss = function(n){
  round(prod(n)/sum(n),0) 
}

power = chars$direct_zeros %>% rename(ctrl = `Treatment Description.x`,
                                      trt = `Treatment Description.y`) %>%  unite(comp,trt,ctrl,sep = " vs ") %>% select(comp, ntot,nstud) %>% mutate(num = seq(1,78,1))

comps = power %>% select(comp,ntot) %>% spread(comp,ntot) #Create wide format of unadjusted ns


#Unadjusted
n = power[grep("vs drops$",power$comp),c(1,2,4)] %>%
  mutate(indirect_n =c(
    #drops_sweet_mult vs drops
    eff_ss(c(comps$`drops_phys vs drops_sweet_mult`,comps$`drops_phys vs drops`)) +
      eff_ss(c(comps$`drops_sweet vs drops_sweet_mult`,comps$`drops_sweet vs drops`))+
      eff_ss(c(comps$`drops_acet30 vs drops_sweet`,comps$`drops_acet30 vs drops`)),
    
    #drops_phys vs drops
    eff_ss(c(comps$`drops_phys vs drops_sweet_mult`,comps$`drops_sweet_mult vs drops`)) +
      eff_ss(c(comps$`drops_sweet vs drops_phys`,comps$`drops_sweet vs drops`))+
      eff_ss(c(comps$`drops_ebm_mult vs drops_sweet_mult`,comps$`drops_ebm_mult vs drops`)),
    
    
    #drops_sweet vs drops
    eff_ss(c(comps$`drops_sweet vs drops_phys`,comps$`drops_phys vs drops`)) +
      eff_ss(c(comps$`drops_sweet vs drops_sweet_mult`,comps$`drops_sweet_mult vs drops`))+
      eff_ss(c(comps$`drops_acet30 vs drops_sweet`,comps$`drops_acet30 vs drops`)),
    
    #placebo vs drops
    0,
    
    #drops_acet30 vs drops
    eff_ss(c(comps$`drops_acet30 vs drops_sweet`,comps$`drops_sweet vs drops`)),
    
    #drops_acet60 vs drops
    0,
    
    #drops_ebm_mult vs drops
    eff_ss(c(comps$`drops_ebm_mult vs drops_sweet_mult`,comps$`drops_sweet_mult vs drops`)) +
      eff_ss(c(comps$`drops_ebm_mult vs drops_phys`,comps$`drops_phys vs drops`)),
    
    
    #drops vs N20 sweet
    eff_ss(c(comps$`drops_N2O_sweet vs drops_sweet`,comps$`drops_sweet vs drops`)),
    
    #drops vs drops WFDRI
    0,
    
    
    #drops vs sweet
    0,
    
    #drops vs sweet_rep
    eff_ss(c(comps$`sweet_rep vs placebo`,comps$`placebo vs drops`)),
    
    #drops vs sweet_sing
    eff_ss(c(comps$`sweet_sing vs placebo`,comps$`placebo vs drops`)))
  )





power_table = n %>% mutate(tot_eff = ntot + indirect_n,
                                                                  sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                                  power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2)) %>% 
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,num)

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

outcome_names = power[grep("vs drops$",power$comp),c(1)] %>% separate(comp,c("t2","t1"), sep = " vs ") %>% select(t2)

pa_reac_pt_res = data.frame(outcome = outcome_names$t2,
                             median = 0,
                             low = 0,
                             high = 0,
                             sd = 0)
for(i in seq_along(outcome_names$t2)){
  temp = summary(relative.effect(pa_reac_data$sa1$results, t1 = "drops",t2 = outcome_names$t2[i])) 
  
  pa_reac_pt_res$median[i] = temp$summaries$quantiles[1,3]
  pa_reac_pt_res$low[i] = temp$summaries$quantiles[1,1]
  pa_reac_pt_res$high[i] = temp$summaries$quantiles[1,5]
  pa_reac_pt_res$sd[i] = temp$summaries$statistics[1,2]
  
  
}


power_table = power_table %>% mutate(post_sd = round(pa_reac_pt_res$sd,2),
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
pa_reac_forest_data = pa_reac_pt_res


pa_reac_forest_data = bind_cols(power_table,pa_reac_forest_data) %>% select(-one_of(c("sample_frac","tot_eff","power","num","outcome","post_sd"))) %>%
  mutate(comp = c("Drops + sweet taste mult","Drops + phys","Drops + sweet taste", "Placebo",
                  "Drops + ebm mult","Drops + acetaminophen 30sec","Drops + acetaminophen 60sec",
                  "Drops + N2O + sweet taste","Drops + WFDRI","Sweet taste alone","Repeated sweet taste",
                  "Sweet taste + singing")) %>% arrange(median) 

pa_reac_plot_data = pa_reac_forest_data %>% select(median, low, high) %>% rename(mean = median,
                                                                                   lower = low,
                                                                                   upper = high)
pa_reac_table_data = pa_reac_forest_data %>%  mutate(mean_cri = paste(round(median,2)," (",round(low,2)," to ",round(high,2),")",sep="")) %>% select(-median,-low,-high,-sd)


pa_reac_table_data$gelman_n = round(pa_reac_table_data$gelman_n,0)
pa_reac_table_data$gelman_m = round(pa_reac_table_data$gelman_m,2)
pa_reac_table_data$gelman_power = round(pa_reac_table_data$gelman_power,2)

pa_reac_table_data = pa_reac_table_data %>%  mutate(sig = c("yes","yes","no","yes","no","yes",rep("no",6)),
                                                    gelman_m = ifelse(sig == "yes",gelman_m,"NA"))  %>% select(comp,ntot,indirect_n,gelman_n,gelman_power,gelman_m,mean_cri)



# pdf("./figs/final/pain scales reactivity/nma_pa_reac_power_forest.pdf", onefile = FALSE, width = 12, height = 5)

png("./figs/final/pain scales reactivity/nma_pa_reac_power_forest.png",res = 300, width = 4000, height = 1700)

forestplot(rbind(c("Comparison","Direct N","Effective 
indirect N
                   ","Heterogeneity
adjusted N
                   ","Power","Exaggeration 
factor
                   ","Mean Difference (95% CrI)"),pa_reac_table_data),
           c(NA,as.numeric(as.character(pa_reac_plot_data$mean))),
           c(NA,as.numeric(as.character(pa_reac_plot_data$lower))),
           c(NA,as.numeric(as.character(pa_reac_plot_data$upper))),
           graph.pos = 7, graphwidth = unit(50,'mm'),
           is.summary = c(TRUE,rep(FALSE,12)),
           hrzl_lines = gpar(col="#444444"),
           align = c("l",rep("c",5),"l"),
           boxsize = 0.4,
           colgap = unit(3,"mm"),
           lineheight = unit(7,"mm"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 11, family = "calibri")),
           col = fpColors(box = "mediumpurple", line = "midnightblue"),
           xlab = "Compared to drops alone",
           title = "PIPP Reactivity (Intervention vs Anesthetic eye drops alone)"
)
dev.off()


save(pa_reac_data, file = "./cache/pa_reac.rda")


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
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# actual timepoint coded as 1 if in immediate (1-2 mins) otherwise 0
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



pa_recov_data = NULL

pa_recov_data$pa$gemtc = prep_gem(pa_recov)

pa_recov_data$pa$gemtc$data = pa_recov_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(pa_recov[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))



pa_recov_data$pa$network = mtc.network(data.re =pa_recov_data$pa$gemtc$data[,1:4],
                                      studies =pa_recov_data$pa$gemtc$data[,c(1,5:10)])

pa_recov_data$pa$results = mtc.model(pa_recov_data$pa$network, type = "consistency",
                                    linearModel = "random",likelihood = "normal",
                                    link = "identity")
pa_recov_data$pa$results = mtc.run(pa_recov_data$pa$results)

summary(pa_recov_data$pa$results)

pa_recov_sucra = sucra(pa_recov_data$pa$results, direction = -1)
pa_recov_sucra = as.data.frame(pa_recov_sucra) %>% rownames_to_column("treatment")

pa_recov_data$pa$anohe = mtc.anohe(pa_recov_data$pa$network)
# plot(summary(pa_recov_data$pa$anohe))

forest(relative.effect.table(pa_recov_data$pa$results),"drops")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- Informative priors on sigma, huge improvement
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data$sa1$results = mtc.model(pa_recov_data$pa$network, type = "consistency",
                                      linearModel = "random",likelihood = "normal",
                                      hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.0625), 
                                      link = "identity")

pa_recov_data$sa1$results = mtc.run(pa_recov_data$sa1$results)



forest(relative.effect.table(pa_recov_data$sa1$results),"drops")
summary(pa_recov_data$sa1$results)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2 
# ----- Meta-regression on cross-over design
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
sa2_regressor = list(coefficient = "shared",
                     variable = "design",
                     control = "drops")

pa_recov_data$sa2$results = mtc.model(pa_recov_data$pa$network, type = "regression",
                                     linearModel = "random",likelihood = "normal",
                                     regressor = sa2_regressor,
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.625),
                                     link = "identity")

pa_recov_data$sa2$results = mtc.run(pa_recov_data$sa2$results)

summary(pa_recov_data$sa2$results)
forest(relative.effect.table(pa_recov_data$sa2$results),"drops")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 
# ----- Meta-regression on timing
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


sa3_regressor = list(coefficient = "shared",
                     variable = "actual_timepoint",
                     control = "drops")



pa_recov_data$sa3$results = mtc.model(pa_recov_data$pa$network, type = "regression",
                                     linearModel = "random",likelihood = "normal",
                                     regressor = sa3_regressor,
                                     link = "identity")

pa_recov_data$sa3$results = mtc.run(pa_recov_data$sa3$results)

summary(pa_recov_data$sa3$results)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3 
# ----- Meta-regression on imputed mean (currently none)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# sa3_regressor = list(coefficient = "shared",
#                      variable = "imputed_mean",
#                      control = "drops")
# 
# 
# 
# pa_recov_data$sa3$results = mtc.model(pa_recov_data$pa$network, type = "regression",
#                                      linearModel = "random",likelihood = "normal",
#                                      regressor = sa3_regressor,
#                                      link = "identity")
# 
# pa_recov_data$sa3$results = mtc.run(pa_recov_data$sa3$results)
# 
# summary(pa_recov_data$sa3$results)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4
# ----- Meta-regression on scaled scores (currently none)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# pa_recov_data$sa4$gemtc = prep_gem(pa_recov %>% filter(imputed_mean == "no"))
# 
# pa_recov_data$sa4$network = mtc.network(data.re =pa_recov_data$sa4$gemtc$data)
# 
# 
# pa_recov_data$sa4$results = mtc.model(pa_recov_data$sa4$network, type = "consistency",
#                                      linearModel = "random",likelihood = "normal",
#                                      link = "identity")
# 
# pa_recov_data$sa4$results = mtc.run(pa_recov_data$sa4$results)
# 
# 
# summary(pa_recov_data$sa4$results)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 5
# ----- Meta-regression on pub type (currently only journal articles)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# pa_recov_data$sa5$gemtc = prep_gem(pa_recov %>% left_join(rop_data_study[c("studlab","pub_type")],by = "studlab") %>% filter(pub_type == "journal"))
# 
# pa_recov_data$sa5$network = mtc.network(data.re =pa_recov_data$sa5$gemtc$data)
# 
# 
# pa_recov_data$sa5$results = mtc.model(pa_recov_data$sa5$network, type = "consistency",
#                                      linearModel = "random",likelihood = "normal",
#                                      link = "identity")
# 
# pa_recov_data$sa5$results = mtc.run(pa_recov_data$sa5$results)
# 
# 
# summary(pa_recov_data$sa5$results)




#=========================================
# Sensitivity 6 - Meta-regression on high overall risk of bias
#=========================================

sa6_regressor = list(coefficient = "shared",
                     variable = "oa_rob_sub",
                     control = "drops")



pa_recov_data$sa6$results = mtc.model(pa_recov_data$pa$network, type = "regression",
                                     linearModel = "random",likelihood = "normal",
                                     hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.625),
                                     regressor = sa6_regressor,
                                     link = "identity")

pa_recov_data$sa6$results = mtc.run(pa_recov_data$sa6$results)

summary(pa_recov_data$sa6$results)


#=========================================
# Sensitivity 8 - Meta-regression on control arm risk
#=========================================

pa_recov_data$sa7 = prep_wb(pa_recov)

pa_recov_data$sa7$meta_cr = pa_recov_data$sa7$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

(pa_recov_data$sa7$list = nma_winbugs_datalist(pa_recov_data$sa7$meta_cr,pa_recov_data$sa7$treatments, contrast = FALSE))

pa_recov_data$sa7$list$mx = as.vector(pa_recov_data$sa7$meta_cr %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


pa_recov_data$sa7$bugs = jags(pa_recov_data$sa7$list,NULL,params_mr,model.file = "./jags_models/re_normal_armdata_meta_jags.txt",
                             n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10)



pa_recov_data$sa7$bugs = nma_outputs(model = pa_recov_data$sa7$bugs$BUGSoutput,pa_recov_data$sa7$treatments)

pa_recov_data$sa7$bugs$B = pa_recov_data$sa7$bugs$bugs[grep("^B$", rownames(pa_recov_data$sa7$bugs$bugs)),]

pa_recov_data$sa7$meta_cr %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                                         y3_diff = y_3 - y_1) %>% gather(diff,value,y2_diff:y3_diff) %>% gather(trt,num, t_2:t_3) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(pa_recov_data$sa7$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)



#==========================================================================
# Power analysis
# Assumptions:
# 1 = I2 is 50%
# 2 = I2 is 70%
# Effect to detect = 2 points (1 MID)
# Assumed SD of PIPP = used Dhaliwhal (largest study)

#==========================================================================

pdf("pa_recov.pdf")
momlinc_netgraph(pa_recov_int,recov_chars$int_char,2)
dev.off()



recov_power = recov_chars$direct_zeros %>% rename(ctrl = `Treatment Description.x`,
                                                  trt = `Treatment Description.y`) %>%  unite(comp,trt,ctrl,sep = " vs ") %>% select(comp, ntot,nstud) %>% mutate(num = seq(1,55,1))

recov_comps = recov_power %>% select(comp,ntot) %>% spread(comp,ntot) #Create wide format of unadjusted ns




#Unadjusted

recov_n = recov_power[grep("vs drops$",recov_power$comp),c(1,2,4)] %>%
  mutate(indirect_n =c(
  
    #drops_phys vs drops
    eff_ss(c(recov_comps$`drops_sweet vs drops_phys`,recov_comps$`drops_sweet_mult vs drops`)) +
      eff_ss(c(recov_comps$`drops_ebm_mult vs drops_phys`,recov_comps$`drops_ebm_mult vs drops`)),
    
    #drops_sweet_mult vs drops
    eff_ss(c(recov_comps$`drops_ebm_mult vs drops_sweet_mult`,recov_comps$`drops_ebm_mult vs drops`)),
    
    #drops_sweet vs drops
    eff_ss(c(recov_comps$`drops_acet30 vs drops_sweet`,recov_comps$`drops_acet30 vs drops`)),
    
    #drops_ebm_mult vs drops
    eff_ss(c(recov_comps$`drops_ebm_mult vs drops_sweet_mult`,recov_comps$`drops_sweet_mult vs drops`)),
    
    #drops vs drops acet 30
    eff_ss(c(recov_comps$`drops_acet30 vs drops_sweet`,recov_comps$`drops_sweet vs drops`)),
    
    #drops vs drops acet 60
    eff_ss(c(recov_comps$`drops_morph vs drops_acet60`,recov_comps$`drops_morph vs drops`)),
    
    #drops vs drops_morph
    eff_ss(c(recov_comps$`drops_morph vs drops_acet60`,recov_comps$`drops_acet60 vs drops`)),
    
    #drops vs phys
    0,
    
    #placebo vs drops
    0,
    
    #drops vs sweet
    0
    
  ))

recov_power_table = recov_n %>% mutate(tot_eff = ntot + indirect_n,
                                                                                    sample_frac = ifelse(tot_eff/req_samp >1,">100",round(tot_eff/req_samp*100,2)),
                                                                                    power = round((power.t.test(tot_eff/2,delta = 2, sd = 2.25))$power,2)) %>% 
  
  select(comp,ntot,indirect_n,tot_eff,sample_frac,power,num)



#Pull posterior SDs, means, and 95% CRi
outcome_names = recov_power[grep("vs drops$",recov_power$comp),c(1)] %>% separate(comp,c("t2","t1"), sep = " vs ") %>% select(t2)

pa_recov_pt_res = data.frame(outcome = outcome_names$t2,
                   median = 0,
                   low = 0,
                   high = 0,
                   sd = 0)
for(i in seq_along(outcome_names$t2)){
temp = summary(relative.effect(pa_recov_data$sa1$results, t1 = "drops",t2 = outcome_names$t2[i])) 

pa_recov_pt_res$median[i] = temp$summaries$quantiles[1,3]
pa_recov_pt_res$low[i] = temp$summaries$quantiles[1,1]
pa_recov_pt_res$high[i] = temp$summaries$quantiles[1,5]
pa_recov_pt_res$sd[i] = temp$summaries$statistics[1,2]


}


recov_power_table = recov_power_table %>% mutate(post_sd = round(pa_recov_pt_res$sd,2),
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
pa_recov_forest_data = pa_recov_pt_res


pa_recov_forest_data = bind_cols(recov_power_table,pa_recov_forest_data) %>% select(-one_of(c("sample_frac","tot_eff","power","num","outcome","post_sd"))) %>%
  
  mutate(comp = c("Drops + phys","Drops + sweet taste mult","Drops + sweet taste","Drops + ebm mult","Drops + acetaminophen 30s",
                  "Drops + acetaminophen 60s","Drops + morphine","Phys alone","Placebo","Sweet taste alone")) %>% arrange(median) 

pa_recov_plot_data = pa_recov_forest_data %>% select(median, low, high) %>% rename(mean = median,
                                                                                   lower = low,
                                                                                   upper = high)
pa_recov_table_data = pa_recov_forest_data %>%  mutate(mean_cri = paste(round(median,2)," (",round(low,2)," to ",round(high,2),")",sep="")) %>% select(-median,-low,-high,-sd)

pa_recov_table_data$gelman_n = round(pa_recov_table_data$gelman_n,0)
pa_recov_table_data$gelman_m = round(pa_recov_table_data$gelman_m,2)
pa_recov_table_data$gelman_power = round(pa_recov_table_data$gelman_power,2)

pa_recov_table_data = pa_recov_table_data %>%  mutate(sig = c("yes",rep("no",9)),
                                                      gelman_m = ifelse(sig == "yes",gelman_m,"NA"))  %>% select(comp,ntot,indirect_n,gelman_n,gelman_power,gelman_m,mean_cri)




png("./figs/final/pain scales recovery/nma_pa_recov_power_forest.png",res = 300, width = 4000, height = 1700)
forestplot(rbind(c("Comparison","Direct N","Effective 
indirect N
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

save(pa_recov_data,file = "./cache/pa_recov.rda")



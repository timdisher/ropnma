library(forestplot)
library(personograph)

source("./analyses/final/heart rate/rop_explore_hr.R")
source("./functions/gemtc test/nma_cont_gemtc.R")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
hr_recov_data = NULL

hr_recov_data$pa$gemtc = prep_gem(hr_recov$data)

hr_recov_data$pa$gemtc$data = hr_recov_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(hr_recov$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis ----
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



hr_recov_data$pa$mod =  set_net(hr_recov_data$pa$gemtc$data, type = "fixed")

calc_n(data = hr_recov_data$pa$gemtc$data,input = hr_recov_data$pa$gemtc$input)


gemtc_diag(hr_recov_data$pa$mod$results)

summary(hr_recov_data$pa$mod$results)

plot(hr_recov_data$pa$mod$network)

#Only loop comes from 3arm trial, consistent by default

hr_recov_data$pa$mod$suc


#League table

hr_recov_data$pa$mod$suc
hr_recov_panames = c("EBM \n multisensory + \n TA",
                      "NNS + TA",
                     "Topical \n Anesthetic (TA)",
                     "Sweet taste \n multisensory + \n TA")


recov_basicp = relative.effect(hr_recov_data$pa$mod$results,t1 = c("drops"),preserve.extra = FALSE)
recov_results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0)
recov_order = hr_recov_data$pa$mod$suc %>% mutate(pub_names = hr_recov_panames)


windows()
league_plot(results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0), order = recov_order, textsize = 3.5)

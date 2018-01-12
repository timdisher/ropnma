#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


source("./analyses/final/cry time/rop_explore_cry.R")
source("./functions/gemtc test/nma_cont_gemtc.R")


# #========================================================================================
# 
# 
# Outcome: Cry Reactivity-----
# 
# #========================================================================================

cry_reac_data = NULL

cry_reac_data$pa$gemtc = prep_gem(cry_reac$data)

cry_reac_data$pa$gemtc$data = cry_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(cry_reac$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         design = ifelse(design == "Parallel",0,1),
         pub_type = ifelse(pub_type == "journal",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1),
         scaled_score = ifelse(scaled_score == "no",0,1),
         actual_timepoint = ifelse(actual_timepoint == "5 min post",0,1)) %>% replace_na(list(actual_timepoint = 0))


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
# --------- Vague priors on sigma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



cry_reac_data$pa$mod = set_net(cry_reac_data$pa$gemtc$data, type = "random")

calc_n(data = cry_reac_data$pa$gemtc$data,input = cry_reac_data$pa$gemtc$input)


gemtc_diag(cry_reac_data$pa$mod$results)

summary(cry_reac_data$pa$mod$results)

#No loops
plot(cry_reac_data$pa$mod$network)

cry_reac_data$pa$mod$suc
cry_reac_panames = c("Sweet taste \n multisensory + \n TA","Topical \n Anesthetic (TA)","NNS + TA",
                     "Acetimanophen 60min \n +TA")


recov_basicp = relative.effect(cry_reac_data$pa$mod$results,t1 = c("drops"),preserve.extra = FALSE)
recov_order = cry_reac_data$pa$mod$suc %>% mutate(pub_names = cry_reac_panames)


windows()
league_plot(results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0), order = recov_order, textsize = 3.5)


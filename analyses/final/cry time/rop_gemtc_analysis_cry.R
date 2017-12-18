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



cry_reac_data$pa$mod = set_net(cry_reac_data$pa$gemtc$data)

calc_n(data = cry_reac_data$pa$gemtc$data,input = cry_reac_data$pa$gemtc$input)


gemtc_diag(cry_reac_data$pa$mod$results)

summary(cry_reac_data$pa$mod$results)

#No loops
plot(cry_reac_data$pa$mod$network)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- Informative priors on sigma, huge improvement
# Large diff = 5 points = 1/25= 0.04
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

cry_reac_data$sa1$results = mtc.model(cry_reac_data$pa$network, type = "consistency",
                                      linearModel = "random",likelihood = "normal",
                                      hy.prior = mtc.hy.prior("std.dev","dhnorm",0,0.04), 
                                      link = "identity")

cry_reac_data$sa1$results = mtc.run(cry_reac_data$sa1$results)

cry_reac_sucra = sucra(cry_reac_data$sa1$results, direction = -1)
cry_reac_sucra = as.data.frame(cry_reac_sucra) %>% rownames_to_column("treatment")

summary(cry_reac_data$sa1$results)
forest(relative.effect.table(cry_reac_data$sa1$results),"drops")
# save(cry_reac,file = "./cache/cry_reac.rda")

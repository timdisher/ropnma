#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


source("./analyses/final/spo2/rop_explore_spo2.R")
source("./functions/gemtc test/nma_cont_gemtc.R")

# #========================================================================================
# 
# 
# Outcome: Sp02 Reactivity-----
# 
# #========================================================================================

os_reac_data = NULL

os_reac_data$pa$gemtc = prep_gem(os_reac$data)

os_reac_data$pa$gemtc$data = os_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","oa_rob_sub","design","pub_type")],by = c("study" = "studlab")) %>%
  left_join(os_reac$data[c("studlab","imputed_mean","scaled_score","actual_timepoint")] %>% replace_na(list(imputed_mean = "no")) %>% distinct(),by = c("study" = "studlab")) %>%
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




os_reac_data$pa$mod = set_net(os_reac_data$pa$gemtc$data)

calc_n(data = os_reac_data$pa$gemtc$data,input = os_reac_data$pa$gemtc$input)



gemtc_diag(os_reac_data$pa$mod$results)

summary(os_reac_data$pa$mod$results)


summary(os_reac_data$pa$mod$results)


os_sucra = sucra(os_reac_data$pa$mod$results, direction = 1)

os_sucra = as.data.frame(os_sucra) %>% rownames_to_column("treat") %>% mutate(treat = paste("d.drops.",treat,sep = "")) %>% rename(sucra = os_sucra)


#No closed loops

plot(os_reac_data$pa$mod$network)

forest(relative.effect.table(os_reac_data$pa$results),"drops")




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
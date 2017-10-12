#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

library(coda)
library(R2WinBUGS)
library(tidyverse)
library(netmeta)
library(stargazer)
library(reshape2)
library(forcats)
library(scales)
library(forestplot)
library(gemtc)
library(reshape2)
library(R2WinBUGS)



source("./analyses/final/adverse events/rop_explore_ae.R")
source("./functions/nma_cont.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
source("./functions/nma_utility_functions.R")
# #========================================================================================
# 
# 
# Outcome: Adverse events Reactivity-----
# 
# #========================================================================================


ae_reac$gemtc$data = as.data.frame(ae_reac$data_connect %>% select(studlab, num_events,sample_size,trt_group) %>%
  rename(study = studlab,
         sampleSize = sample_size,
         treatment = trt_group,
         responders = num_events))

ae_reac$pa$network = mtc.network(data.ab =ae_reac$gemtc$data)


plot(ae_reac$pa$network, edge.arrow.size = .4, vertex.size = 30, arrow.mode = 3)

ae_reac$pa$results = mtc.model(ae_reac$pa$network, type = "consistency",
                                     linearModel = "fixed",likelihood = "binom",
                                     link = "logit")

ae_reac$pa$results = mtc.run(ae_reac$pa$results)

forest(relative.effect.table(ae_reac$pa$results),"drops")





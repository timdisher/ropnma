#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

source("./analyses/final/adverse events/rop_explore_ae.R")

source("./functions/gemtc test/nma_cont_gemtc.R")
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

sum(ae_reac$gemtc$data$sampleSize)
length(ae_reac$gemtc$data$sampleSize)

#No loops
plot(ae_reac$pa$network, edge.arrow.size = .4, vertex.size = 30, arrow.mode = 3)

ae_reac$pa$results = mtc.model(ae_reac$pa$network, type = "consistency",
                                     linearModel = "random",likelihood = "binom",
                                     link = "logit")

ae_reac$pa$results = mtc.run(ae_reac$pa$results)

gemtc_diag(ae_reac$pa$results)


summary(ae_reac$pa$results)
ae_sucra = sucra(ae_reac$pa$results, direction = -1)

ae_sucra = as.data.frame(ae_sucra) %>% rownames_to_column("treat") %>% 
  mutate(treat = paste("d.drops.",treat,sep = "")) %>% rename(sucra = ae_sucra)

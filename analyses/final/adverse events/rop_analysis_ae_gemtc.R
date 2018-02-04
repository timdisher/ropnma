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
plot(ae_reac$pa$network)

ae_reac$pa$results = mtc.model(ae_reac$pa$network, type = "consistency",
                                     linearModel = "fixed",likelihood = "binom",
                                     link = "logit")

ae_reac$pa$results = mtc.run(ae_reac$pa$results)

gemtc_diag(ae_reac$pa$results)


summary(ae_reac$pa$results)
ae_sucra = sucra(ae_reac$pa$results, direction = -1)

ae_sucra = as.data.frame(ae_sucra) %>% rownames_to_column("treat") %>% 
  mutate(treat = paste("d.drops.",treat,sep = "")) %>% rename(sucra = ae_sucra) %>% arrange(-sucra)

ae_sucra 
ae_reac_panames = c("Sweet taste + \n TA","Topical \n Anesthetic (TA)", "Acetaminophen 60min + \nTA","No treatment"
                    )


recov_basicp = relative.effect(ae_reac$pa$results,t1 = c("drops"),preserve.extra = FALSE)
recov_order = ae_sucra  %>% mutate(pub_names = ae_reac_panames)


windows()
league_plot(results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0), order = recov_order, textsize = 3.5, exp = TRUE)


# #========================================================================================
# 
# 
# Outcome: Adverse events Regulation-----
# 
# #========================================================================================


ae_recov$gemtc$data = as.data.frame(ae_recov$data %>% select(studlab, num_events,sample_size,trt_group) %>%
                                     rename(study = studlab,
                                            sampleSize = sample_size,
                                            treatment = trt_group,
                                            responders = num_events))


ae_recov$pa$network = mtc.network(data.ab =ae_recov$gemtc$data)

sum(ae_recov$gemtc$data$sampleSize)
length(ae_recov$gemtc$data$sampleSize)

#No loops
plot(ae_recov$pa$network, edge.arrow.size = .4, vertex.size = 30, arrow.mode = 3)

ae_recov$pa$results = mtc.model(ae_recov$pa$network, type = "consistency",
                               linearModel = "fixed",likelihood = "binom",
                               link = "logit")

ae_recov$pa$results = mtc.run(ae_recov$pa$results)

gemtc_diag(ae_recov$pa$results)


summary(ae_recov$pa$results)
recov_ae_sucra = sucra(ae_recov$pa$results, direction = -1)

recov_ae_sucra = as.data.frame(recov_ae_sucra) %>% rownames_to_column("treat") %>% 
  mutate(treat = paste("d.drops.",treat,sep = "")) %>% rename(sucra = recov_ae_sucra) %>% arrange(-sucra)

recov_ae_sucra 
ae_recov_panames = c("No treatment", "Sweet taste + \nTA","Topical \n Anesthetic (TA)"
)


recov_basicp = relative.effect(ae_recov$pa$results,t1 = c("drops"),preserve.extra = FALSE)
recov_order = ae_sucra  %>% mutate(pub_names = ae_recov_panames)


windows()
league_plot(results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0), order = recov_order, textsize = 3.5, exp = TRUE)

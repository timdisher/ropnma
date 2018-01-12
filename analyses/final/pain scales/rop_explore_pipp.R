# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Exploratory analysis and table of characteristics generation for
# retinopathy of prematurity NMA - pain scale outcome.


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Pain reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

source("./analyses/final/rop_import.R")
source("./functions/gemtc test/nma_cont_gemtc.R")





# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------

pa_reac = rop_data_arm %>% filter(outcome == "PIPP" | outcome == "NIPS", timepoint_group == "reactivity") #O'sullivan only study with N-PASS score scaled to PIPP


#excluded trials
(pa_reac_excluded = pa_reac %>% filter(studlab %in% c("Ucar 2014","Olsson 2011")) %>% group_by(studlab) %>% summarize(sample = sum(sample_size),
                                                                                                                     treat = paste(trt_group,collapse=' vs ')) %>% 
  mutate(reason = c("Does not use speculum","No variance info"))
)


#Prepare data

pa_reac = pa_reac %>% filter(!studlab %in% pa_reac_excluded$studlab) #No variance info, Olsson does not use speculum
pa_reac = pa_reac %>% droplevels()


### converts data to correct format to allow for assessment of connectivity

(pa_reac_net = prep_gem(pa_reac))
pa_reac_net = mtc.network(data.re = pa_reac_net$data[,1:4])
plot(pa_reac_net)

#Generate characteristics tables
chars = netmeta_xl_chars(data = pa_reac,outcome = "pa_reac",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph

pa_reac_netgraph = pub_netgraph(chars, nodecolour = c("#7467aa"))
png("./figs/pain scales reactivity/reac_netgraph.png", height = 8.8,width = 10, units = "in", res = 300)
pa_reac_netgraph 
dev.off()

#for pich pres

# pdf("./figs/pain scales reactivity/reac_netgraph_pich.pdf", height = 8.8,width = 10)
# pub_netgraph(chars, nodecolour = c("#7467aa")) + theme(legend.position = "none")
# dev.off()
# 
png("./figs/pain scales reactivity/reac_netgraph_pich2.png", height = 8.8,width = 10, unit = "in", res = 300)
pub_netgraph(chars, nodecolour = c("#7467aa")) + theme(legend.position = "none")
dev.off()



### Placebo response graph
pa_reac_graph =  pairwise(data = pa_reac,treat = trt_group, n= sample_size, mean = mean,sd = std_dev,studlab = studlab) %>%
  mutate(TE = TE*-1) %>% left_join(pa_reac[c("studlab","imputed_mean","scaled_score")] %>% distinct(), by = "studlab") %>%left_join(rop_data_study[c("studlab","design","avg_pma","oa_rob_sub","avg_bw","avg_ga","pub_type")], by = "studlab") %>% 
  mutate(oa_rob_sub = ifelse(oa_rob_sub == "low","low","high"))

pich_graph = pa_reac_graph %>% filter(treat1 == "drops")


treatment_names = c(
  "drops" = "Drops alone",
  "drops_acet30" = "Drops + \n acetaminophen (30s)",
  "drops_acet60" = "Drops + acetaminophen (60s)",
  "drops_ebm_mult" = "Drops + \n EBM multisensory",
  "drops_N2O_sweet" = "Drops + \n N2O + sweet taste",
  "drops_phys" = "Drops + \n physical",
  "drops_sweet" = "Drops + \n sweet taste",
  "drops_sweet_mult" = "Drops + \n sweet multisensory",
  "drops_WFDRI" = "Drops + \n WFDRI",
  "placebo" = "Placebo",
  "sweet" = "Sweet taste \n alone",
  "sweet_rep" = "Repeated sweet \n taste",
  "sweet_sing" = "Sweet + \n singing"
)


### Covariate plots

reg_te_plots = function(data = pa_reac_graph, labeller = treatment_names){

plots = list(
  ctrl = nma_bargraphs(data = data, outcome = "TE", scat = T, x = "mean1", labeller = labeller) + labs(y = "Mean difference (PIPP)", x = "Control arm response") + ggtitle(""),
  
  bw = nma_bargraphs(data = data %>% mutate(avg_bw = ifelse(avg_bw == 0,NA,avg_bw)), labeller = labeller, outcome = "TE", scat = T, x = "avg_bw") + labs(y = "Mean difference (PIPP)", x = "Birthweight") + ggtitle(""),
  
  ga = nma_bargraphs(data = data %>% mutate(avg_bw = ifelse(avg_bw == 0,NA,avg_ga)), labeller = labeller, outcome = "TE", scat = T, x = "avg_ga") + labs(y = "Mean difference (PIPP)", x = "Gestational age") + ggtitle(""),
  
  pma = nma_bargraphs(data = data %>% mutate(avg_bw = ifelse(avg_bw == 0,NA,avg_pma)), labeller = labeller, outcome = "TE", scat = T, x = "avg_pma") + labs(y = "Mean difference (PIPP)", x = "Postmentrual age") + ggtitle(""),
  
  
  # Categorical
  
  im = nma_bargraphs(data = data, outcome = "TE", x = "imputed_mean", scat = T, labeller = labeller) + 
    geom_point(shape = 1,stroke = 1.3, aes(size = seTE, colour = imputed_mean)) + 
    scale_colour_manual(name = "Imputed mean", values = c("darkorchid","lightblue"),labels = c("No","Yes")) +
    labs(y = "Mean difference (PIPP)", x = "Imputed mean") + ggtitle(""),
  
  xover = nma_bargraphs(data = data,outcome = "TE", x = "design", scat = T, labeller = labeller) + 
    geom_point(shape = 1,stroke = 1.3, aes(size = seTE, colour = design)) + 
    scale_colour_manual(name = "Design", values = rev(c("darkorchid","lightblue")),labels = c("Crossover","Parallel")) +
    labs(y = "Mean difference (PIPP)", x = "Design") + ggtitle(""),
  
  scaled = nma_bargraphs(data = data,outcome = "TE", x = "scaled_score", scat = T, labeller = labeller) + 
    geom_point(shape = 1,stroke = 1.3, aes(size = seTE, colour = scaled_score)) + 
    scale_colour_manual(name = "Scaled score", values = c("darkorchid","lightblue"),labels = c("No","Yes")) + labs(y = "Mean difference (PIPP)", x = "Scaled score") + ggtitle(""),
  
  rob = nma_bargraphs(data = data,outcome = "TE", x = "oa_rob_sub", scat = T, labeller = labeller) + 
    geom_point(shape = 1,stroke = 1.3, aes(size = seTE, colour = oa_rob_sub)) + 
    scale_colour_manual(name = "Risk of bias", values = rev(c("darkorchid","lightblue")),labels = c("High","Low")) + labs(y = "Mean difference (PIPP)", x = "Risk of bias") + ggtitle(""),
  
  pubtype = nma_bargraphs(data = data,outcome = "TE", x = "pub_type", scat = T, labeller = labeller) + 
    geom_point(shape = 1,stroke = 1.3, aes(size = seTE, colour = pub_type)) + 
    scale_colour_manual(name = "Publication type", values = c("darkorchid","lightblue"),labels = c("Journal","Poster")) + labs(y = "Mean difference (PIPP)", x = "Publication type") + ggtitle("")
)
plots
}


pa_reac_scats = reg_te_plots()
pdf("./figs/pain scales reactivity/pa_reac_scats.pdf", height = 10, width = 11)
invisible(lapply(pa_reac_scats,print))
dev.off()


pa_reac_scats_pich = reg_te_plots(pich_graph)
pdf("./figs/pain scales reactivity/pa_reac_scats_pich.pdf", height = 4.5, width = 10)
invisible(lapply(pa_reac_scats_pich,print))
dev.off()



# 
# pich_graph

#
#
#
#====== Recovery=====
#
#
#

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------

pa_recov = rop_data_arm %>% filter(outcome == "PIPP", timepoint_group == "recovery") 


### converts data to correct format to allow for assessment of connectivity
pa_recov_contrast = pairwise(data = pa_recov,treat = trt_group, n= sample_size, mean = mean,sd = std_dev,studlab = studlab) 

#- Assess whether network is connected -#
(pa_recov_net = prep_gem(pa_recov))
pa_recov_net = mtc.network(data.re = pa_recov_net$data[,1:4])
plot(pa_recov_net)


#Generate characteristics tables
recov_chars = netmeta_xl_chars(pa_recov,"pa_recov",ref = "drops",treat = "trt_group",location = "./tables/final/char_tables/pain scales recovery") ##create table of characteristics


### Generate network graph
pa_recov_netgraph = pub_netgraph(recov_chars, nodecolour = c("#7467aa"))

png("./figs/pain scales recovery/recov_netgraph.png", height = 8.8,width = 10, units = "in", res = 300)
pa_recov_netgraph
dev.off()


#PICH Network diagram


#pdf("./figs/pain scales recovery/recov_netgraph_pich.pdf", height = 8.8,width = 10)
png("./figs/pain scales recovery/recov_netgraph_pich2.png", height = 8.8,width = 10, units = "in", res = 300)
pub_netgraph(recov_chars, nodecolour = c("#7467aa"))
dev.off()

### Placebo response graph
pa_recov_graph =  pairwise(data = pa_recov,treat = trt_group, n= sample_size, mean = mean,sd = std_dev,studlab = studlab) %>%
  mutate(TE = TE*-1) %>% left_join(pa_recov[c("studlab","imputed_mean","scaled_score")] %>% distinct(), by = "studlab") %>%left_join(rop_data_study[c("studlab","design","avg_pma","oa_rob_sub","avg_bw","avg_ga")], by = "studlab") %>% 
  replace_na(replace = list(avg_bw = 0, imputed_mean = "no",scaled_score = "no",avg_bw = 0)) %>% mutate(oa_rob_sub = ifelse(oa_rob_sub == "low","low","high"))



treatment_names_recov = c(
  "drops" = "Drops alone",
  "drops_acet30" = "Drops + \n acetaminophen (30s)",
  "drops_acet60" = "Drops + acetaminophen (60s)",
  "drops_ebm_mult" = "Drops + \n EBM multisensory",
  "drops_phys" = "Drops + \n physical",
  "drops_sweet" = "Drops + \n sweet taste",
  "drops_sweet_mult" = "Drops + \n sweet multisensory",
  "placebo" = "Placebo",
  "sweet" = "Sweet taste \n alone",
  "drops_morph" = "Drops + \n morphine",
  "phys" = "Physical \n intervention alone"
  
)

### meta-regression plots

pa_recov_scats = reg_te_plots(pa_recov_graph, labeller = treatment_names_recov)

pdf("./figs/pain scales recovery/pa_recov_scats.pdf", height = 10, width = 11)
invisible(lapply(pa_recov_scats,print))
dev.off()



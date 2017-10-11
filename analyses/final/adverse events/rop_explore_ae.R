# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Exploratory analysis and table of characteristics generation for
# retinopathy of prematurity NMA - Adverse event scale outcome.


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Adverse event reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

source("./analyses/final/rop_import.R")
library(netmeta)
library(forcats)
library(stargazer)
library(stringr)
library(grid)




# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Adverse event reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
ae_reac = NULL

ae_reac$overview = rop_data_arm %>% filter(outcome %in% c("any_ae","02 desat <85%","tachy >180bpm","brady < 100bpm","02 desat > 10%","02 desat >10% on PIPP points",
                                                         "brady and desat","Apnea in 24h","02 < 88% in 24h",
                                                         "desat < 80%"), timepoint_group == "reactivity") #No variability data for saunders, mehta is presence/absence



# As per cochrane, we can analyse change cores and raw scores together. 

(ae_reac$data = ae_reac$overview %>% filter(outcome %in% c("any_ae"))) 

### converts data to correct format to allow for assessment of connectivity
ae_reac$contrast = pairwise(data = ae_reac$data,treat = trt_group, n= sample_size, 
                             event = num_events,studlab = studlab, sm = "OR") 



#- Assess whether network is connected -#
(ae_reac$netconnect = netconnection(treat1,treat2,data = ae_reac$contrast) ) ###Two networks, but one of them is a single study (Dilli)

ae_reac$data_connect = ae_reac$data %>% filter(studlab != "Dilli 2014")

ae_reac_excluded = ae_reac$overview %>% filter(studlab == "Dilli 2014") %>%
  select(studlab,outcome,sample_size) %>% group_by(studlab,outcome) %>% summarise(n = sum(sample_size)) %>% ungroup() %>%
  mutate(reason = c(rep("not part of connected network",4)))

ae_reac$contrast_nodili = pairwise(data = ae_reac$data_connect,treat = trt_group, n= sample_size, 
                            event = num_events,studlab = studlab, sm = "OR") 

(ae_reac$netconnect_nodilli = netconnection(treat1,treat2,data = ae_reac$contrast_nodili) ) ###Two networks, but one of them is a single study (Dilli)



(ae_reac$int = netmeta(TE,seTE,treat1,treat2,studlab,data = ae_reac$contrast_nodili, sm = "OR", comb.random = TRUE)) 


#Generate characteristics tables
(ae_reac$chars = netmeta_xl_chars(ae_reac$data_connect,"ae",ref = "drops",treat = "trt_group", cont = FALSE, event = "num_events")) ##create table of characteristics


### Generate network graph
momlinc_netgraph(ae_reac$int,ae_reac$chars$int_char,2)


### Placebo response graph
(ae_reac$graph$data = ae_reac$data %>% left_join(rop_data_study, by = c("studlab")))

(ae_reac$graph$pr = plac_resp_graph(ae_reac$data, ref = "drops"))

(ae_reac$graph$pr_spec = plac_resp_graph(ae_reac$graph$data,y = "num_events", ref = "drops", facet = TRUE, fv = "speculum"))


# All pairwise meta-analyses

ae_reac$graph$pw = all_pairwise(ae_reac$chars$direct,ae_reac$contrast_nodili, outcome = "ae", 
                                cont = FALSE,sm = "OR",location = "./figs/final/ae_reac/pairwise ma")


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Adverse event recovery----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
ae_recov = NULL

(ae_recov$overview = rop_data_arm %>% filter(outcome %in% c("any_ae","02 desat <85%","tachy >180bpm","brady < 100bpm","02 desat > 10%","02 desat >10% on PIPP points",
                                                          "brady and desat","Apnea in 24h","02 < 88% in 24h",
                                                          "desat < 80%"), timepoint_group == "recovery")) #No variability data for saunders, mehta is presence/absence



# As per cochrane, we can analyse change cores and raw scores together. 

(ae_recov$data = ae_recov$overview %>% filter(outcome %in% c("any_ae"))) 

### converts data to correct format to allow for assessment of connectivity 
ae_recov$contrast = pairwise(data = ae_recov$data,treat = trt_group, n= sample_size, 
                            event = num_events,studlab = studlab, sm = "OR") 



#No connected networks, below code kept for future analyses.


# #- Assess whether network is connected -#
# (ae_recov$netconnect = netconnection(treat1,treat2,data = ae_recov$contrast) ) ###Two networks, but one of them is a single study (Dilli)
# 
# ae_recov$data_connect = ae_recov$data %>% filter(studlab != "Dilli 2014")
# 
# ae_recov$contrast_nodili = pairwise(data = ae_recov$data_connect,treat = trt_group, n= sample_size, 
#                                    event = num_events,studlab = studlab, sm = "OR") 
# 
# (ae_recov$netconnect_nodilli = netconnection(treat1,treat2,data = ae_recov$contrast_nodili) ) ###Two networks, but one of them is a single study (Dilli)
# 
# 
# 
# (ae_recov$int = netmeta(TE,seTE,treat1,treat2,studlab,data = ae_recov$contrast_nodili, sm = "OR", comb.random = TRUE)) 
# 
# 
# #Generate characteristics tables
# (ae_recov$chars = netmeta_xl_chars(ae_recov$data_connect,"ae",ref = "drops",treat = "trt_group", cont = FALSE, event = "num_events")) ##create table of characteristics
# 
# 
# ### Generate network graph
# momlinc_netgraph(ae_recov$int,ae_recov$chars$int_char,2)
# 
# 
# ### Placebo response graph
# (ae_reac$graph$data = ae_reac$data %>% left_join(rop_data_study, by = c("studlab")))
# 
# (ae_reac$graph$pr = plac_resp_graph(ae_reac$data, ref = "drops"))
# 
# (ae_reac$graph$pr_spec = plac_resp_graph(ae_reac$graph$data,y = "num_events", ref = "drops", facet = TRUE, fv = "speculum"))
# 
# 
# # All pairwise meta-analyses
# 
# ae_reac$graph$pw = all_pairwise(ae_reac$chars$direct,ae_reac$contrast_nodili, outcome = "ae", 
#                                 cont = FALSE,sm = "OR",location = "./figs/final/ae_reac/pairwise ma")
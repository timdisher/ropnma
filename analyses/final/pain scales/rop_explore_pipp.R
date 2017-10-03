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
library(netmeta)
library(forcats)
library(stargazer)
library(stringr)
library(grid)

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------

pa_reac = rop_data_arm %>% filter(outcome == "PIPP" | outcome == "NIPS", timepoint_group == "reactivity") #O'sullivan only study with N-PASS score scaled to PIPP

pa_reac = pa_reac %>% filter(studlab != "Ucar 2014") #No variance info
pa_reac$trt_group = fct_infreq(pa_reac$trt_group) %>% droplevels()

### converts data to correct format to allow for assessment of connectivity
pa_reac_contrast = pairwise(data = pa_reac,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 


#- Assess whether network is connected -#
pa_reac_netconnect = netconnection(treat1,treat2,data = pa_reac_contrast)
pa_reac_int = netmeta(TE,seTE,treat1,treat2,studlab,data = pa_reac_contrast, sm = "MD", comb.random = TRUE) ###required to drawn netgraph

#Generate characteristics tables
chars = netmeta_xl_chars(pa_reac,"pa_reac",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph
momlinc_netgraph(pa_reac_int,chars$int_char,2)


### Placebo response graph
pa_reac_graph = pa_reac %>% left_join(rop_data_study, by = c("studlab"))

plac_resp_graph(pa_reac, ref = "drops")

plac_resp_graph(pa_reac_graph, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

(pa_reac_pairwise = all_pairwise(chars$direct,pa_reac_contrast, outcome = "PIPP reactivity", location = "./figs/final/pain scales reactivity/pairwise ma", sm = "MD"))


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Pain reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------

pa_recov = rop_data_arm %>% filter(outcome == "PIPP", timepoint_group == "recovery") 

### converts data to correct format to allow for assessment of connectivity
pa_recov_contrast = pairwise(data = pa_recov,treat = trt_group, n= sample_size, mean = mean,sd = std_dev,studlab = studlab) 

#- Assess whether network is connected -#
(pa_recov_netconnect = netconnection(treat1,treat2,data = pa_recov_contrast))

(pa_recov_int = netmeta(TE,seTE,treat1,treat2,studlab,data = pa_recov_contrast, sm = "MD")) ###required to drawn netgraph

#Generate characteristics tables
recov_chars = netmeta_xl_chars(pa_recov,"pa_recov",ref = "drops",treat = "trt_group",location = "./tables/final/char_tables/pain scales recovery") ##create table of characteristics


### Generate network graph
momlinc_netgraph(pa_recov_int,recov_chars$int_char,2)

### Placebo response graph
pa_recov_graph = pa_recov %>% left_join(rop_data_study, by = c("studlab"))

plac_resp_graph(pa_recov, ref = "drops")

plac_resp_graph(pa_recov_graph, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

(pa_recov_pairwise = all_pairwise(recov_chars$direct,sm = "MD",pa_recov_contrast, outcome = "PIPP recovery", location = "./figs/final/pain scales recovery/pairwise ma"))


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




# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Pain reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
hr_reac = NULL

hr_reac$overview = rop_data_arm %>% filter(grepl("hr",outcome), timepoint_group == "reactivity") 

# Can use MD for change from baseline, and raw heart rate. Mehta will have to be dropped because there is no variance and it's a crossover. Max heart rate is a different outcome,
# can't combine it plus change scores because it would require standardizing which cochrane does not recommend.

hr_reac$data = hr_reac$overview %>% filter(studlab %in% c("Dhaliwal 2010","Olsson 2011","Xin 2016"), outcome != "hr max") #Drops Mehta and connects network. Will hgave to report taplak seperately.

### converts data to correct format to allow for assessment of connectivity
hr_reac$contrast = pairwise(data = hr_reac$data,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 



#- Assess whether network is connected -#
hr_reac$netconnect = netconnection(treat1,treat2,data = hr_reac$contrast) #two sub-networks, no way to connect. One network 
# is entirely supported by a single multi arm study. Decision is to report this seperate.

hr_reac$int = netmeta(TE,seTE,treat1,treat2,studlab,data = hr_reac$contrast, sm = "MD", comb.random = TRUE) ###required to drawn netgraph

#Generate characteristics tables
hr_reac$chars = netmeta_xl_chars(hr_reac$data,"hr",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph
hr_reac$graph$network = momlinc_netgraph(hr_reac$int,hr_reac$chars$int_char,2)


### Placebo response graph
hr_reac$graph$data = hr_reac$data %>% left_join(rop_data_study, by = c("studlab"))

hr_reac$graph$pr = plac_resp_graph(hr_reac$data, ref = "drops")

hr_reac$graph$pr_spec = plac_resp_graph(hr_reac$graph$data, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

hr_reac$graph$pw = all_pairwise(hr_reac$chars$direct,hr_reac$contrast, outcome = "heart rate", sm = "SMD",location = "./figs/final/hr_reac/pairwise ma")


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----Pain recovery----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
hr_recov = NULL

hr_recov$overview = rop_data_arm %>% filter(grepl("hr",outcome), timepoint_group == "recovery") 

# Challenge here as always is mix of outcomes. Hr max, change from baseline, and raw heart rate. Decision will be to use
# standardized mean difference.Will have to report Mehta seperately because they provide no variance AND are a crossover.

hr_recov$data = hr_recov$overview %>% filter(studlab != "Mehta 2005", outcome != "hr max") #Drops Mehta and connects network. Will have to report taplak seperately.

### converts data to correct format to allow for assessment of connectivity
hr_recov$contrast = pairwise(data = hr_recov$data,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 



#- Assess whether network is connected -#
hr_recov$netconnect = netconnection(treat1,treat2,data = hr_recov$contrast) #two sub-networks, no way to connect. One network 
# is entirely supported by a single multi arm study. Decision is to report this seperate.

hr_recov$int = netmeta(TE,seTE,treat1,treat2,studlab,data = hr_recov$contrast, sm = "MD", comb.random = TRUE) ###required to drawn netgraph

#Generate characteristics tables
hr_recov$chars = netmeta_xl_chars(hr_recov$data,"hr",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph
hr_recov$graph$network = momlinc_netgraph(hr_recov$int,hr_recov$chars$int_char,2)


### Placebo response graph
hr_recov$graph$data = hr_recov$data %>% left_join(rop_data_study, by = c("studlab"))

hr_recov$graph$pr = plac_resp_graph(hr_recov$data, ref = "drops")

hr_recov$graph$pr_spec = plac_resp_graph(hr_recov$graph$data, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

hr_recov$graph$pw = all_pairwise(hr_recov$chars$direct,hr_recov$contrast, outcome = "heart rate", sm = "MD",location = "./figs/final/hr_recov/pairwise ma")



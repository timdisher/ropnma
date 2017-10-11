# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Exploratory analysis and table of characteristics generation for
# retinopathy of prematurity NMA - oxygen saturation scale outcome.


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----oxygen saturation reactivity----
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
#                    ----oxygen saturation reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
os_reac = NULL

os_reac$overview = rop_data_arm %>% filter(grepl("02",outcome), timepoint_group == "reactivity") 


# As per cochrane, we can analyse change cores and raw scores together. 

os_reac$data = os_reac$overview %>% filter(outcome %in% c("02 decrease","02 sat"), !(studlab %in% c("Olsson 2011"))) %>% 
  mutate(mean = ifelse(outcome == "02 decrease",mean*-1,mean))#Drops Mehta and connects network. Will hgave to report taplak seperately.


os_reac_excluded = os_reac$overview %>% filter(!(outcome %in% c("02 decrease","02 sat")) | studlab == "Olsson 2011") %>%
  select(studlab,outcome,sample_size) %>% group_by(studlab,outcome) %>% summarise(n = sum(sample_size)) %>% ungroup() %>%
  mutate(reason = c("not an outcome of interest",
                    "analyzed with adverse events",
                    "analyzed with adverse events",
                    "not an outcome of interest",
                    "no variance",
                    "no speculum"))



### converts data to correct format to allow for assessment of connectivity
os_reac$contrast = pairwise(data = os_reac$data,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 



#- Assess whether network is connected -#
os_reac$netconnect = netconnection(treat1,treat2,data = os_reac$contrast) 

os_reac$int = netmeta(TE,seTE,treat1,treat2,studlab,data = os_reac$contrast, sm = "MD", comb.random = TRUE) ###required to drawn netgraph

#Generate characteristics tables
os_reac$chars = netmeta_xl_chars(os_reac$data,"sp02",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph
momlinc_netgraph(os_reac$int,os_reac$chars$int_char,2)


### Placebo response graph
os_reac$graph$data = os_reac$data %>% left_join(rop_data_study, by = c("studlab"))

os_reac$graph$pr = plac_resp_graph(os_reac$data, ref = "drops")

os_reac$graph$pr_spec = plac_resp_graph(os_reac$graph$data, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

os_reac$graph$pw = all_pairwise(os_reac$chars$direct,os_reac$contrast, outcome = "sp02", sm = "MD",location = "./figs/final/os_reac/pairwise ma")


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----oxygen saturation recovery----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
os_recov = NULL

os_recov$overview = rop_data_arm %>% filter(grepl("02",outcome), timepoint_group == "recovery") 

# As per cochrane, we can analyse change cores and raw scores together. 

os_recov$data = os_recov$overview %>% filter(outcome %in% c("02 sat change","02 sat"), studlab != "Mehta 2005") 


os_recov_excluded = os_recov$overview %>% filter(!(outcome %in% c("02 decrease","02 sat")) | studlab == "Mehta 2005") %>%
  select(studlab,outcome,sample_size) %>% group_by(studlab,outcome) %>% summarise(n = sum(sample_size)) %>% ungroup() %>%
  mutate(reason = c("analyzed with adverse events",
                    "analyzed with adverse events",
                    "not an outcome of interest",
                    "no variance"))

### converts data to correct format to allow for assessment of connectivity
os_recov$contrast = pairwise(data = os_recov$data,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 



#- Assess whether network is connected -#
os_recov$netconnect = netconnection(treat1,treat2,data = os_recov$contrast) 

os_recov$int = netmeta(TE,seTE,treat1,treat2,studlab,data = os_recov$contrast, sm = "MD", comb.random = TRUE) ###required to drawn netgraph

#Generate characteristics tables
os_recov$chars = netmeta_xl_chars(os_recov$data,"sp02",ref = "drops",treat = "trt_group") ##create table of characteristics


### Generate network graph
momlinc_netgraph(os_recov$int,os_recov$chars$int_char,2)


### Placebo response graph
os_recov$graph$data = os_recov$data %>% left_join(rop_data_study, by = c("studlab"))

os_recov$graph$pr = plac_resp_graph(os_recov$data, ref = "drops")

os_recov$graph$pr_spec = plac_resp_graph(os_recov$graph$data, ref = "drops", facet = TRUE, fv = "speculum")


# All pairwise meta-analyses

os_recov$graph$pw = all_pairwise(os_recov$chars$direct,os_recov$contrast, outcome = "sp02", sm = "MD",location = "./figs/final/os_recov/pairwise ma")


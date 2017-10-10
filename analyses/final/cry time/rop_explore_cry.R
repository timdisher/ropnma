# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Exploratory analysis and table of characteristics generation for
# retinopathy of prematurity NMA - crying time scale outcome.


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----crying time reactivity----
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
#                    ----crying time reactivity----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
cry_reac = NULL

cry_reac$overview = rop_data_arm %>% filter(grepl("cry",outcome), timepoint_group == "reactivity") #No variability data for saunders, mehta is presence/absence


# As per cochrane, we can analyse change cores and raw scores together. 

cry_reac$data = cry_reac$overview %>% filter(outcome %in% c("crying time"),!(trt_group %in% c("drops_diet_1hr","drops_diet_2hr")),!(studlab %in% c("Olsson 2011"))) %>%
  mutate(std_dev = ifelse(is.na(std_dev),se*sqrt(sample_size),std_dev)) #diet not in connected network

### converts data to correct format to allow for assessment of connectivity
cry_reac$contrast = pairwise(data = cry_reac$data,treat = trt_group, n= sample_size, 
                            mean = mean,sd = std_dev,studlab = studlab) 


cry_excluded = data.frame(study = c("Mehta 2005",
                                    "Strube 2010",
                                    "Olsson 2011"),
                          reason = c("Outcome is presence/absence",
                                     "Diet is not in connected network",
                                     "No speculum used"))
#- Assess whether network is connected -#
(cry_reac$netconnect = netconnection(treat1,treat2,data = cry_reac$contrast) )

(cry_reac$int = netmeta(TE,seTE,treat1,treat2,studlab,data = cry_reac$contrast, sm = "MD", comb.random = TRUE)) ###required to drawn netgraph

#Generate characteristics tables
(cry_reac$chars = netmeta_xl_chars(cry_reac$data,"sp02",ref = "drops",treat = "trt_group")) ##create table of characteristics


### Generate network graph
momlinc_netgraph(cry_reac$int,cry_reac$chars$int_char,2)


### Placebo response graph
(cry_reac$graph$data = cry_reac$data %>% left_join(rop_data_study, by = c("studlab")))

(cry_reac$graph$pr = plac_resp_graph(cry_reac$data, ref = "drops"))

(cry_reac$graph$pr_spec = plac_resp_graph(cry_reac$graph$data, ref = "drops", facet = TRUE, fv = "speculum"))


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                    ----crying time recovery----
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Generate tables and graphs----------------------------------------------------------------------------------------------------------------------------------------------
cry_recov = NULL

cry_recov$overview = rop_data_arm %>% filter(grepl("cry",outcome), timepoint_group == "recovery") #No variability data for saunders, mehta is presence/absence

View(cry_recov$overview)

#Only mehta
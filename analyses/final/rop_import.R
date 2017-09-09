library(tidyverse)
source("./functions/netmetaxl_tables_function.R")

(rop_data_study = read_csv("./data/rop_studylevel_aug.csv")[-1,]) #Exclude variable description row

(rop_data_arm = read_csv("./data/rop_armlevel_aug.csv")[-1,])

(codebook = read_csv("./data/rop_var_types.csv"))

(int_codes = read_csv("./data/rop_intcodes.csv"))


###As this is the final data analysis, we replace all NDs and SNs with NAs.

rop_data_study = replace(rop_data_study,rop_data_study == "ND" | rop_data_study == "SN", NA)

rop_data_arm = replace(rop_data_arm,rop_data_arm == "ND" | rop_data_arm == "SN", NA)


###Recode int classes 
(rop_data_arm$trt_group = int_codes[match(rop_data_arm$treatment,int_codes[["intervention_name"]]),2][[1]])

###Use codebook to transform variables to correct type

codebook_study = codebook %>% filter(sheet == "study_level")

rop_data_study = lookup_type(rop_data_study,codebook_study)


codebook_arm = codebook %>% filter(sheet == "arm_level")

rop_data_arm = lookup_type(rop_data_arm,codebook_arm)



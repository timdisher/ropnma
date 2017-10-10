prep_gem = function(data){
  data$trt_group = fct_infreq(data$trt_group) %>% droplevels(data)#reorders factor from most to least common
  
  data = droplevels(data) # Drops factor levels that don't exist (otherwise they are carried over)
  long_gemtc(data = data)
  
}


#============================================================================================
#
#
#             ---- Arrange long data into winbugs format -----
#
#
#============================================================================================

#data = pa_reac, a data frame limited to the outcome of interest



long_gemtc = function(data = pa_reac,studlab = "studlab", trt = "trt_group",mean = "mean",sd = "std_dev",sample_size = "samplesize",studylevel = rop_data_study){
  
  (input = data %>% select(studlab,trt_group,mean,std_dev,sample_size,p_value) %>% 
     rename(y = mean,
            n = sample_size,
            sd = std_dev) %>%  
     arrange(studlab,trt_group) %>% ##Arranges treatments by factor level
     group_by(studlab) %>%
     mutate(arm = row_number()) %>% rename(treatment = trt_group)
  )
  
 
data = input %>% group_by(studlab) %>% mutate(y = y - first(y)) %>% mutate(y = ifelse(arm == 1,NA,y)) %>%
    mutate(std.err = ifelse(arm == 1, round(sd/sqrt(n),2),
                            se_md(first(sd),sd,first(n),n))) %>% left_join(studylevel %>% select(studlab,design), by = studlab) %>%
  mutate(std.err = ifelse(design == "Crossover" & arm != 1,se_paired(y,p_value,first(n)),std.err)) %>%
  rename(diff = y,
         study = studlab) %>% select(study,treatment,diff,std.err) %>% as.data.frame()##Did not take square of V as in netmeta xl, revisit if issues
    
  ###Output a list
  
list(input = input,data = data)
}









#===================================================================================================













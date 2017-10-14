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
  mutate(std.err = ifelse(design == "Crossover" & arm != 1 & !is.na(p_value),se_paired(y,p_value,first(n)),std.err)) %>% mutate(narm = max(arm)) %>%
  mutate(std.err = ifelse(narm == 2 & arm == 1, NA, std.err)) %>%
  rename(diff = y,
         study = studlab) %>% select(study,treatment,diff,std.err) %>% as.data.frame()
    
  ###Output a list
  
list(input = input,data = data)
}

library(ggnetwork)
library(sna)


pub_netgraph = function(chars, data, nodecolour = "mediumorchid2", layout = "circle"){
  
  chars_net = chars$int_char %>% select(Treatment, `Number of Patients`) %>% rename(sample = `Number of Patients`)
  
  df = chars$direct %>% rename(weight = `# Studies`) %>% select(Treatment, weight) %>% separate(Treatment,c("from","to"), sep = " vs ")
  df = graph.data.frame(df)
  V(df)$size = as.numeric(chars_net[match(V(df)$name,chars_net$Treatment),2][[1]])
  
  legend = data.frame(treat = V(df)$name, letter = LETTERS[as.factor(V(df)$name)]) %>% mutate(legend = paste(letter,": ",treat,sep = "")) %>% select(legend) %>%
    arrange(legend)
  
  network = ggnetwork(df, layout = layout)
  
  ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(size = weight), color = "grey50", alpha = 0.8) +
    geom_nodes(aes(size = size,colour = vertex.names)) + 
    theme_blank() + 
    geom_edgelabel_repel(aes(label = weight), color = "grey25") +
    geom_nodetext(aes(label = LETTERS[vertex.names]))+
    scale_size_area("sample", max_size = 32,guide = FALSE) + scale_colour_manual(name = "Treatment", values = rep(nodecolour,length(network$x)),
                                                                                 labels = legend[[1]]) + theme(legend.position = "bottom")
  
  
}




#===================================================================================================


#Calculate SUCRA----
# Adapted from Dias WinBUGS codes

sucra = function(results,direction){
  rank = rank.probability(results, preferredDirection = direction)
  rank = rank[,1:length(rank[1,])]
  cumeffectiveness = rank
  SUCRA = rank[,1]
  for(k in seq_along(rank[,1])){
    for(h in seq_along(rank[1,])){
      cumeffectiveness[k,h] = sum(rank[k,1:h])
    }
  }
  
  for(i in seq_along(cumeffectiveness[,1])){
    SUCRA[i] = sum(cumeffectiveness[i, 1:(length(cumeffectiveness[,1]) - 1)])/(length(cumeffectiveness[,1]) - 1)
  }
  SUCRA
}






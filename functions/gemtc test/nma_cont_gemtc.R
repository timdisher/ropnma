library(tidyverse)
library(forcats)
library(stargazer)
library(stringr)
library(grid)
library(reshape2)
library(netmeta)
library(gemtc)
library(gridExtra)
library(R2jags)
library(scales)
library(forestplot)

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

(wide= input %>%
    recast(studlab ~ variable + arm, id.var = c(studlab,"arm"))
)
  
list(input = input,data = data, wide = wide)
}

library(ggnetwork)
library(sna)
library(igraph)
library(intergraph)

pub_netgraph = function(chars,nodecolour = "mediumorchid2", layout = "circle"){
  
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
                                                                                 labels = legend[[1]]) + theme(legend.position = "bottom",
                                                                                                               legend.text = element_text(size = 14))
  

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


#=========== =
# Interrfacing with JAGS/WB #====

prep_wb = function(data, smd = FALSE,convert_contrast = TRUE){
  data$trt_group = fct_infreq(data$trt_group) %>% droplevels(data)#reorders factor from most to least common
  
  data = droplevels(data) # Drops factor levels that don't exist (otherwise they are carried over)
  
  
  if(convert_contrast == TRUE){
    if(smd == TRUE) long_wb_smd(data, hedge = TRUE) else long_wb(data = data)}else long_arm_wb(data = data)
  
}


long_wb = function(data = pa_reac,studlab = "studlab", trt = "trt_group",mean = "mean",sd = "std_dev",sample_size = "samplesize"){
  
  (input = data %>% select(studlab,trt_group,mean,std_dev,sample_size) %>% 
     rename(y = mean,
            n = sample_size,
            sd = std_dev) %>%  
     arrange(studlab,trt_group) %>% ##Arranges treatments by factor level
     group_by(studlab) %>%
     mutate(arm = row_number(),
            t = as.numeric(trt_group),
            na = n())
  )
  
  treatments = input %>% ungroup() %>% select(trt_group,t) %>% distinct() %>% arrange(t) %>%
    rename(description = trt_group)
  
  (wide= input %>% select(-trt_group) %>%
      recast(studlab ~ variable + arm, id.var = c(studlab,"arm")) %>% select(studlab:na_1) %>% rename(na = na_1)
  )
  
  
  if("y_4" %in% colnames(wide)){
    
    wide_con = wide %>% mutate(y_2 = y_2 - y_1,
                               y_3 = y_3 - y_1,
                               y_4 = y_4 - y_1,
                               se_2 = se_md(sd_1,sd_2,n_1,n_2),
                               se_3 = se_md(sd_1,sd_3,n_1,n_3),
                               se_4 = se_md(sd_1,sd_4,n_1,n_4),
                               V = ifelse(na > 2,(sd_1/sqrt(n_1))^2, NA)) %>% arrange(na) #arrange by number of arms (see WiBUGS code for why)
    
    winbugs_ready = wide_con %>% select(t_1:t_4,y_2:y_4,se_2:se_4,V,na)} else if("y_3" %in% colnames(wide)){
      
      
      wide_con = wide %>% mutate(y_2 = y_2 - y_1,
                                 y_3 = y_3 - y_1,
                                 se_2 = se_md(sd_1,sd_2,n_1,n_2),
                                 se_3 = se_md(sd_1,sd_3,n_1,n_3),
                                 V = ifelse(na > 2,(sd_1/sqrt(n_1))^2, NA)) %>% arrange(na) #arrange by number of arms (see WiBUGS code for why)
      
      winbugs_ready = wide_con %>% select(t_1:t_3,y_2:y_3,se_2:se_3,V,na)} else {
        
        wide_con = wide %>% mutate(y_2 = y_2 - y_1,
                                   se_2 = se_md(sd_1,sd_2,n_1,n_2)) %>% arrange(na) #arrange by number of arms (see WiBUGS code for why)
        
        winbugs_ready = wide_con %>% select(t_1:t_2,y_2,se_2,na)
        
      }
  
  #####Output a list
  
  out = list(input = input,treatments = treatments,arm_wide = wide,wide = wide_con,wb = winbugs_ready)
}


#####Long to wb for smd data. V is calculated differently here. This iteration was found
# http://methods.cochrane.org/sites/methods.cochrane.org.cmi/files/public/uploads/S8-L%20Problems%20introduced%20by%20multi-arm%20trials%20-%20full%20network%20meta-analysis.pdf

long_wb_smd= function(data = pa_reac,studlab = "studlab", trt = "trt_group",mean = "mean",sd = "std_dev",sample_size = "samplesize", hedge = FALSE){
  
  (input = data %>% select(studlab,trt_group,mean,std_dev,sample_size) %>% 
     rename(y = mean,
            n = sample_size,
            sd = std_dev) %>%  
     arrange(studlab,trt_group) %>% ##Arranges treatments by factor level
     group_by(studlab) %>%
     mutate(arm = row_number(),
            t = as.numeric(trt_group),
            na = n())
  )
  
  treatments = input %>% ungroup() %>% select(trt_group,t) %>% distinct() %>% arrange(t) %>%
    rename(description = trt_group)
  
  (wide= input %>% select(-trt_group) %>%
      recast(studlab ~ variable + arm, id.var = c(studlab,"arm")) %>% select(studlab:na_1) %>% rename(na = na_1)
  )
  
  if("y_4" %in% colnames(wide)){
    
    (wide = wide %>% mutate(y_2 = smd(y_2,y_1,sd_2,sd_1,n_2,n_1),
                            y_3 = smd(y_3,y_1,sd_3,sd_1,n_3,n_1),
                            y_4 = smd(y_4,y_1,sd_4,sd_1,n_4,n_1)))
    
    (wide = wide %>%  mutate(se_2 = se_smd(y_2,sd_2,sd_1,n_2,n_1),
                             se_3 = se_smd(y_3,sd_3,sd_1,n_3,n_1),
                             se_4 = se_smd(y_4,sd_4,sd_1,n_4,n_1),
                             V = ifelse(na > 2,1/n_1,NA)) %>% arrange(na)) } else if("y_3" %in% colnames(wide)){
                               
                               
                               (wide = wide %>% mutate(y_2 = smd(y_2,y_1,sd_2,sd_1,n_2,n_1),
                                                       y_3 = smd(y_3,y_1,sd_3,sd_1,n_3,n_1)))
                               
                               (wide = wide %>%  mutate(se_2 = se_smd(y_2,sd_2,sd_1,n_2,n_1),
                                                        se_3 = se_smd(y_3,sd_3,sd_1,n_3,n_1),
                                                        V = ifelse(na > 2,1/n_1,NA)) %>% arrange(na)) } else {
                                                          
                                                          
                                                          (wide = wide %>% mutate(y_2 = smd(y_2,y_1,sd_2,sd_1,n_2,n_1)))
                                                          
                                                          (wide = wide %>%  mutate(se_2 = se_smd(y_2,sd_2,sd_1,n_2,n_1),
                                                                                   V = ifelse(na > 2,1/n_1,NA)) %>% arrange(na))}
  
  
  #arrange by number of arms (see WiBUGS code for why)
  if(hedge == FALSE){
    
    wide %>% select(studlab,matches("t_"),matches("y_"),matches("se_"),V,na) %>% select(-y_1)
  } else{
    
    
    out =  wide %>% select(studlab,matches("t_"),matches("y_"),matches("se_"),V,na) %>% select(-y_1)
    list(wblist = out, treatments = treatments)
    
  }
}


#================
# Create data list
#===============

nma_winbugs_datalist = function(data,treatments,contrast = TRUE){
  
  
  #Find column indexes
  t_loc = grep("t_|t. != studlab",colnames(data),fixed = F)
  y_loc = grep("y_|y.",colnames(data),fixed = F)
  se_loc = grep("se_|se.",colnames(data),fixed = F)
  na_loc = grep("na",colnames(data),fixed = F)
  
  if(contrast == TRUE){
    v_loc = grep("V",colnames(data),fixed = F)
    
    
    t = as.matrix(data[,t_loc])
    
    # number of treatments
    nt = length(treatments$description)
    
    y = as.matrix(cbind(rep(NA, length(t[,1])), data[,y_loc]))
    
    se = as.matrix(cbind(rep(NA, length(t[,1])), data[,se_loc]))
    
    na = as.vector(data[,na_loc])
    
    V = as.vector(data[,v_loc])
    
    ns2 = length(subset(na, na==2))
    ns3 = length(subset(na, na==3))
    ns4 = length(subset(na, na==4))
    
    if(ns4 >0)  data = list(nt=nt, ns2=ns2, ns3=ns3, ns4=ns4,t=t, y=y, se=se, na=na, V=V) else if(ns3 >0) data = list(nt=nt, ns2=ns2, ns3=ns3,t=t, y=y, se=se, na=na, V=V) else data = list(nt=nt, ns2=ns2,t=t, y=y, se=se, na=na)
    
    
    
    data
  } else {t = as.matrix(data[,t_loc])
  
  nt = length(treatments$description)
  
  y = as.matrix(data[,y_loc])
  
  se = as.matrix(data[,se_loc])
  
  na = as.vector(data[,na_loc])
  
  ns = length(t[,1])
  
  data = list(nt=nt,t=t, y=y, se=se, na=na, ns = ns)
  
  data}
  
  
}

#============================================================
# Create outputs - For use when doing meta-regression
#============================================================
nma_outputs = function(model, treatments,FE = FALSE){
  #======================================
  # Output all pairwise comparisons
  #======================================
  
  # summarize mean differences - mean, 2.5%, median, 97.5%
  md = (as.matrix(model$summary[grep("meandif", rownames(model$summary), fixed=F), c(1,3,5,7)]))
  
  # summarize the probability better - mean, 2.5%, median, 97.5%
  better = (as.matrix(model$summary[grep("better", rownames(model$summary), fixed=F), c(1,2)]))
  
  
  # conver comparison columns into a matrix
  b <- gsub("\\meandif|\\[|\\]", "",rownames(md))
  b_clean <- matrix(as.numeric(unlist(strsplit(b, ","))), ncol=2, byrow=T)
  comp <- paste(treatments$description[b_clean[,2]], ' vs.', treatments$description[b_clean[,1]])
  
  md_med <- round(md[,"50%"],2)
  
  md_ci <- paste(paste(round(md[,'2.5%'],2), round(md[,'97.5%'],2), sep =" to "), sep="")
  
  hr_med_ci = paste(round(md[,"50%"],2) , paste("(", paste(round(md[,'2.5%'],2), round(md[,'97.5%'],2), sep =" to "), ")", sep=""))
  
  better_mean_sd = paste(round(better[,"mean"],2) , paste("(", paste(round(better[,'sd'],2), ")", sep=""), sep=""))
  
  better_mean = round(better[,"mean"],2)
  
  better_sd = round(better[,'sd'],2)
  
  
  comp = cbind(comp, md_med, md_ci, better_mean, better_sd); colnames(comp) = c('Comparison (Trt A vs. Trt B)', 'Mean Difference of Trt A vs. Trt B', "95% CrI of Mean Difference", 'Probability of Trt A better than Trt B', "SD of Probability Better")
  rownames(comp) <- NULL
  
  
  #======================================
  #SUCRA, Rankings, and Probability Best Plots
  #======================================
  
  
  # sucra summary
  sucra = model$summary[grep("SUCRA", rownames(model$summary), fixed=F),] #mean SUCRA, remove ",1" to get median values
  
  # ranks summary
  ranks = model$summary[grep("rk", rownames(model$summary), fixed=F),] #mean rank, remove ",1" to get median values
  
  # prob best summary
  probs = model$summary[grep("best", rownames(model$summary), fixed=F),] #mean probability best, remove ",1" to get median values
  
  # mean of ranks (2.5% of rank to 97.5% of rank)
  rk_mean_ci = paste(round(ranks[,"50%"],2) , paste("(", paste(round(ranks[,'2.5%'],2), round(ranks[,'97.5%'],2), sep =" to "), ")", sep=""))
  
  orders = c('Mean SUCRA',"Mean Probability of Best Treatment", "Mean Rank (95% CrI)")
  rankings = cbind(100*sucra[,1],     # mean of SUCRA, in %
                   100*probs[,1],     # mean of prob best, in % 
                   as.vector(rk_mean_ci))     # mean of ranks (2.5% of rank to 97.5% of rank)
  colnames(rankings) = orders
  rownames(rankings) = treatments$description
  
  rr = rankings
  
  
  #Probability summary
  probs <- model$mean$prob; rownames(probs) <- treatments$description; colnames(probs) <- seq(1,length(treatments$description))
  
  dat <- melt(cbind(probs)); colnames(dat) <- c('Treatments', 'Rankings', 'Probability of Best Treatment')
  
  
  rankogram = ggplot(dat,aes(x = Rankings, y = `Probability of Best Treatment`,fill = Treatments)) + 
    geom_bar(position = "fill",stat = "identity") + 
    scale_y_continuous(labels = percent_format())  + scale_x_continuous(breaks=seq(1,length(treatments$description)))
  
  
  
  if(FE == TRUE){  
    results = list(model = model,comp = comp,rr = rr,rankogram = rankogram,bugs = model$summary, dic = model$DIC)
    
  }else {sd = model$summary[grep("^sd$",row.names(model$summary),fixed = F),]
  
  results = list(model = model,sd = sd,comp = comp,rr = rr,rankogram = rankogram,bugs = model$summary, dic = model$DIC)}
  
  
  results
}




#These functions are commonly used calculations for systematic reviews======



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#
#             ---- Impute mean from Median -----
# Adapted from AHRQ report on handling continuous outcomes in quantitative synthesis
# https://ahrq-ehc-application.s3.amazonaws.com/media/pdf/continuous-outcomes-quantitative-synthesis_methods.pdf
#
# IQR based method adapted from Wiebe et al "A systematic review identifies a lack of standardization in methods
# for handling missing variance data" Journal of Clinical Epidemilogy 50(2006)
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#Min and max----
mean_med = function(min,median,max){
  (min + 2*median + max)/4
}


sd_med = function(min,median,max){
  sqrt(1/12*(((min-2*median)^2/4)+(max - min)^2))
}

#IQR----
sd_med_iqr = function(upper, lower){
  (upper-lower)/1.35
}



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#
#             ---- Calculate the standard error of a mean difference -----
#
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

se_md = function(sd1,sd2,n1,n2){
  
  se = sqrt((sd1^2*(n1-1) + sd2^2*(n2-1))/(n1+n2-2))*sqrt(1/n1+1/n2)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#
#             ---- Calculate standard error from a paired t-test -----
#
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

se_paired = function(y_2,p_value,n_1){
  se = abs(y_2/qt((p_value/2),n_1-1))
  se}


#============================================================================================
#
#
#             ---- Calculate standardized mean difference  -----
# Source: https://cran.r-project.org/web/packages/compute.es/compute.es.pdf
# Per recommendation from Sofia Dias - will need to adapt this to use pooled sd across all arms
#============================================================================================

smd = function(yE,yC,sdE,sdC,nE,nC){
  n = nE+nC
  
  
  sd = sqrt(((nE-1)*(sdE^2)+(nC-1)*(sdC^2))/(n-2))
  
  md = yE-yC
  
  (md/sd) * (1 - (3/(4*n-9)))
  
}

se_smd = function(g,sdE,sdC,nE,nC){
  n = nE+nC
  
  
  var = n/(nE*nC) + (g^2/(2*(n-3.94)))
  
  sqrt(var)
}


sd_smd = function(yE,yC,sdE,sdC,nE,nC){
  sqrt(((nE-1)*(sdE^2)+(nC-1)*(sdC^2))/(n-2))
}



# Table functions =================



#-----------------------------------------------------------------------------------------------------------------------------

#- Do before function -#

#-----------------------------------------------------------------------------------------------------------------------------


# pa_reac = data_arm %>% filter(outcome == "PIPP" | outcome == "NIPS" | outcome == "N-PASS", timepoint.group == "reactivity") ##choose outcome
# trt_names = c("drops",
#             "plac",
#             "sweet.drops",
#             "nns.drops",
#             "sweet.nns.drops",
#             "wfdri.drops",
#             "paracetamol.drops",
#             "NIDCAP",
#             "NO.sweet.drops",
#             "sweet",
#             "EBM.drops",
#             "ss.drops",
#             "nns.ebm.drops",
#             "sweet.sing",
#             "sweet.rep")## provide a list of names
# 
# ### Example entry
# netmeta_xl_chars(filter(pa_reac,outcome != "N-PASS"),"test")


#-----------------------------------------------------------------------------------------------------------------------------
###- Do in function ###-
##requires tidyverse and netmeta
#-----------------------------------------------------------------------------------------------------------------------------

netmeta_xl_chars = function(data,outcome,ref,treat = "treatment",location = getwd(),cont = TRUE,event = NULL){
  
  if(cont == TRUE){  
    contrast = pairwise(data = data,treat = data[[treat]], n= sample_size, mean = mean,sd = std_dev,studlab = studlab)} else{
      
      contrast = pairwise(data = data,treat = data[[treat]], n= sample_size, event = data[[event]],studlab = studlab)
    }
  
  
  
  studies = (data %>% select(studlab) %>% distinct() %>% count())[[1]] ## count the number of studies
  trts = (data %>% select(treat) %>% distinct() %>% count())[[1]] ## count the number of treatments
  totn = data %>% select(sample_size) %>% sum() # count the number of patients
  poss_pw = choose(trts,2) # total possible pairwise comparisons
  
  
  #-- Create a numbered list of treatments with ref as #1--#
  data[[treat]] = fct_infreq(data[[treat]])
  data = droplevels(data)
  names= (data %>% select(treat) %>% distinct())[[treat]]
  
  names = tibble("Treatment Number" = as.numeric(names),
                 "Treatment Description" = names) %>% arrange(`Treatment Number`)
  
  
  #- Assign treatment numbers based on tibble above -#
  contrast = names %>% rename(treat2_asnum = `Treatment Number`,
                              treat2 = `Treatment Description`) %>% right_join(contrast, by = "treat2")
  contrast = names %>% rename(treat1_asnum = `Treatment Number`,
                              treat1 = `Treatment Description`) %>% right_join(contrast, by = "treat1")
  
  
  
  
  #- Counts number of comparisons with direct evidence -#
  direct = tibble(trt1 = combn(trts,2)[1,],
                  trt2 = combn(trts,2)[2,],
                  nstud = 0,
                  ntot = 0)
  
  for(i in 1:length(direct$trt1)){
    temp = contrast %>% filter(treat1_asnum == direct$trt1[i] & treat2_asnum == direct$trt2[i] | 
                                 treat1_asnum == direct$trt2[i] & treat2_asnum == direct$trt1[i]) %>% 
      summarise(n_stud = n(),
                sample = sum(n1,n2))
    
    direct$nstud[i] = temp[[1]]
    direct$ntot[i] = temp[[2]]
  }
  tot_direct = direct %>% filter(nstud > 0) %>% count()
  
  
  ###calculate number of arms
  arm_table = table(contrast$studlab)
  narms = tibble(test = names(arm_table),
                 vals = arm_table) %>% 
    mutate(twoarm = ifelse(vals <3,1,0),
           threearm = ifelse(vals == 3,1,0),
           fourarm = ifelse(vals == 6,1,0))
  
  twoarm = sum(narms$twoarm)
  multarm = studies - twoarm
  
  #-----------------------------------
  # Network characteristics table
  #-----------------------------------
  net_char <- tibble("Characteristic" = c("Number of Interventions",
                                          "Number of Studies",
                                          "Total Number of Patients in Network",
                                          "Total Possible Pairwise Comparisons",
                                          "Total Number Pairwise Comparisons With Direct Data",
                                          "Number of Two-arm Studies",
                                          "Number of Multi-arms Studies"),
                     "Number" = c(trts,
                                  studies,
                                  totn,
                                  poss_pw,
                                  tot_direct[[1]],
                                  twoarm,
                                  multarm
                     ))
  
  
  #-----------------------------------
  # Intervention characteristics table
  #-----------------------------------,
  style = "asq"
  int_nums = data %>% group_by_(treat) %>% summarise("Number of Comparisons" = n(),
                                                     "Number of Patients" = sum(sample_size)) %>% rename_(Treatment = treat)
  
  
  int_char <- tibble("Treatment" = names[[2]]) %>% left_join(int_nums, by = "Treatment")
  
  
  
  #-----------------------------------
  # Direct comparison characteristics table
  #-----------------------------------
  direct_comp_char <- direct %>% filter(nstud > 0) %>% left_join(names, by = c("trt1" = "Treatment Number")) %>%
    left_join(names, by = c("trt2" = "Treatment Number")) %>% rename(treatment_1 = `Treatment Description.x`,
                                                                     treatment_2 = `Treatment Description.y`) %>%
    mutate(Treatment = paste(treatment_2,"vs",treatment_1)) %>%
    select(Treatment,nstud,ntot) %>% rename("# Studies" = nstud,
                                            "# Patients" = ntot)
  
  direct = direct %>% left_join(names, by = c("trt1" = "Treatment Number")) %>%
    left_join(names, by = c("trt2" = "Treatment Number"))
  
  out = list(net_char = net_char, int_char = int_char, direct = direct_comp_char,nma_data = contrast,direct_zeros = direct)
  
  out
}




#-----------------------------------------------------------------------------------------------------
# This function uses code from netmeta to ensure that treatment codes are organized in ascending order
# As written, it takes results from the pairwise function.Ultimately will need to think of how this works for winbugs
#-----------------------------------------------------------------------------------------------------
treat_reorg = function(){
  wo = as.character(pa_reac_contrast$treat1) > as.character(pa_reac_contrast$treat2)
  
  pa_reac_contrast$TE[wo] = -pa_reac_contrast$TE[wo]
  ttreat1 = as.character(pa_reac_contrast$treat1)
  pa_reac_contrast$treat1[wo] = pa_reac_contrast$treat2[wo]
  pa_reac_contrast$treat2[wo] = ttreat1[wo]
}


#-----------------------------------------------------------------------------------------------------
# This function formats data type according to a data dictionary
#-----------------------------------------------------------------------------------------------------
# codes = subset(factor_codes,sheet = "study")
# test = data_study_orig

lookup_type = function(data,codes,name,date_format = "%m/%d/%Y"){
  
  for(i in seq_along(data)){
    if(is.na(codes[match(colnames(data[i]),codes[[1]]),4])){
      paste(colnames(data[i]),"is missing from the coding diary")} else
        
        if(codes[match(colnames(data[i]),codes[[1]]),4] == "factor"){
          data[[i]] = as.factor(data[[i]])
        } else
          
          if(codes[match(colnames(data[i]),codes[[1]]),4] == "dbl"){
            data[[i]] = as.numeric(data[[i]])
          } else
            
            if(codes[match(colnames(data[i]),codes[[1]]),4] == "date"){
              data[[i]] = as.Date(data[[i]], paste(date_format))
            } else              
              
              data[[i]] = data[[i]]
            
  }
  data
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# This function converts numbers to words, and can be used for manuscript writing
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r
numbers2words <- function(x){
  ## Function by John Fox found here: 
  ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  ## Tweaks by AJH to add commas and "and"
  helper <- function(x){
    
    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) as.vector(ones[digits])
    else if (nDigits == 2)
      if (x <= 19) as.vector(teens[digits[1]])
    else trim(paste(tens[digits[2]],
                    Recall(as.numeric(digits[1]))))
    else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and", 
                                      Recall(makeNumber(digits[2:1]))))
    else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(Recall(makeNumber(digits[
        nDigits:(3*nSuffix + 1)])),
        suffixes[nSuffix],"," ,
        Recall(makeNumber(digits[(3*nSuffix):1]))))
    }
  }
  trim <- function(text){
    #Tidy leading/trailing whitespace, space before comma
    text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
    #Clear any trailing " and"
    text=gsub(" and$","",text)
    #Clear any trailing comma
    gsub("\ *,$","",text)
  }  
  makeNumber <- function(...) as.numeric(paste(..., collapse=""))     
  #Disable scientific notation
  opts <- options(scipen=100) 
  on.exit(options(opts)) 
  ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine") 
  names(ones) <- 0:9 
  teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
             "sixteen", " seventeen", "eighteen", "nineteen")
  names(teens) <- 0:9 
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
            "ninety") 
  names(tens) <- 2:9 
  x <- round(x)
  suffixes <- c("thousand", "million", "billion", "trillion")     
  if (length(x) > 1) return(trim(sapply(x, helper)))
  helper(x)
}

###=========================NMA BARGRAPHS====
nma_bargraphs = function(data = pa_reac_graph, outcome = "avg_pma", labeller = "label_value", scat = F, x = "studlab"){
 
  
  mean_line = data %>% select_("treat1",outcome) %>% group_by(treat1) %>% summarize_all(funs(mean(.,na.rm = TRUE))) %>% rename_(mean = outcome)
  
  if(scat == F){
  data %>% ggplot(aes_string(x = x, y = outcome)) + 
    geom_hline(data = mean_line, aes(yintercept = mean,linetype = "Mean"), size = 1, colour = c("#00CD66")) +
    geom_col(fill = "darkorchid",width = 0.8, position = "dodge") +
    theme_bw() +
    labs(title = outcome,
         x = "Study ID",
         y = outcome)  + 
    facet_grid(treat1~treat2, scales = "free_x", space = "free", switch = "both", labeller = labeller(treat1 = labeller, treat2 = labeller)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 13),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.spacing = unit(0.3,"lines"),
          strip.placement = "outside",
          strip.text.x = element_text(angle = 90),
          strip.text.y = element_text(angle = 180),
          strip.background = element_blank(),
          panel.background = element_blank(),
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())  +
    scale_y_continuous(expand = c(0,0)) + scale_linetype_manual(name = "", values = 1)} else{
  
  
  
  data %>% ggplot(aes_string(x = x, y = outcome)) +
    geom_point(colour = "darkorchid", shape = 1, stroke = 1.3, aes(size = seTE)) + geom_smooth(method = "lm", colour = c("#00CD66"), se = FALSE) +
    theme_bw() +
    labs(title = outcome,
         x = "Study ID",
         y = outcome)  + 
    facet_grid(treat1~treat2, space = "free", switch = "both", labeller = labeller(treat1 = labeller, treat2 = labeller)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 13),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.spacing = unit(0.8,"lines"),
          strip.placement = "outside",
          strip.text.x = element_text(angle = 90),
          strip.text.y = element_text(angle = 180),
          strip.background = element_blank(),
          panel.background = element_blank(),
          
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())  + scale_size(guide = FALSE)
  }
  
  
}

#Diagnostic plots
gemtc_diag = function(results = pa_reac_data$pa$results){
  windows(record = TRUE)
  plot(results)
  gelman.plot(results)
  gelman.diag(results)
  
}


#Create a league table

league = function(results = pa_reac_basicp, order = pa_reac_sucra){
  
  t = results %>% select(order$treat)
  
  df = matrix(nrow = length(t), ncol = length(t))
  
  order = order %>% separate(treat, into = c("a","treat"), sep = "d.drops.") %>% select(-a)
  
  for(i in 1:length(t)){
    
  cons = t[,i] - t
  
  temp = purrr::map(cons, ~quantile(.,c(0.025,0.5,0.975))) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column()
  
  df[,i] = paste(round(temp[,3],2)," (",round(temp[,2],2)," to ",round(temp[,4],2),")",sep = "")
  
  df[i,i] = order[i,1]
  }
  df[upper.tri(df)] = NA
  df
}


# League table plot =========================

league_plot = function(results = pa_reac_basicp, #Basic parameters (d's)
                       order = pa_reac_sucra,  #What is to be used for ordering (SUCRA in default)
                       pub = TRUE, #If using different names for publication
                       colour = "#d2cee3", #Light purple for NAs
                       textsize = 3, #Size for text
                       exp = FALSE, #If d's are on log scale
                       statsig = 0 #Threshold for statsig (change to 1 if ORs)
){
  
  t = results %>% select(order$treat)
  
  
  #Populate required data frame
  df = data.frame(lwr=NA, med=NA, upr=NA, trt1=NA, trt2=NA)
  
  
  #Fill will all comparisons
  for(i in 1:length(t)){
    
    cons =  t[,i] - t 
    
    temp = purrr::map(cons, ~quantile(.,c(0.025,0.5,0.975))) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column()
    
    trt1 = colnames(t[i])
    temp = temp %>% filter(rowname != trt1)
    
    df = df %>% add_row(lwr = temp[,2],med = temp[,3],upr = temp[,4], trt1 = trt1, trt2 = temp$rowname)
  }
  
  
  
  
  if(pub == TRUE){
    df["trt1"] = order[match(df[["trt1"]],order[["treat"]]),3]
    df["trt2"] = order[match(df[["trt2"]],order[["treat"]]),3]
    
    df = df %>% mutate(trt1 = factor(df$trt1, levels = order[["pub_names"]]),
                       trt2 = factor(df$trt2, levels = order[["pub_names"]]))
    if(exp == TRUE){ 
      df$val = ifelse(is.na(df$lwr),NA,paste(round(exp(df$med), 2), "\n (", round(exp(df$lwr), 2), " to ", round(exp(df$upr), 2),")", sep = ""))
      
    } else{
      
      df$val = ifelse(is.na(df$lwr),NA,paste(round(df$med, 2), "\n (", round(df$lwr, 2), " to ", round(df$upr, 2),")", sep = ""))
    }
    
    df = df %>% mutate(val = ifelse(as.numeric(trt1) > as.numeric(trt2), NA, val)) %>% mutate(trt1 = factor(df$trt1, levels = order[["pub_names"]]),
                                                                                              trt2 = factor(df$trt2, levels = rev(order[["pub_names"]])),
                                                                                              sig = ifelse(is.na(val),1,ifelse(lwr < statsig & upr < statsig | lwr > statsig & upr > statsig,2,0))) #Identify stat sig results
    
    ext = data_frame(lwr=NA, med=NA, upr=NA, trt1=order[["pub_names"]], trt2=order[["pub_names"]], val=order[["pub_names"]], sig = NA)
    
    df = rbind(df,ext) %>% slice(-1)
    
    
    
  }  else if(pub == FALSE){
    
    
    df = df %>% mutate(trt1 = factor(df$trt1, levels = order$treat),
                       trt2 = factor(df$trt2, levels = order$treat))
    
    
    if(exp == TRUE){ 
      df$val = ifelse(is.na(df$lwr),NA,paste(round(exp(df$med), 2), "\n (", round(exp(df$lwr), 2), " to ", round(exp(df$upr), 2),")", sep = ""))
      
    } else{
      
      df$val = ifelse(is.na(df$lwr),NA,paste(round(df$med, 2), "\n (", round(df$lwr, 2), " to ", round(df$upr, 2),")", sep = ""))
    }
    
    
    
    df = df %>% mutate(val = ifelse(as.numeric(trt1) > as.numeric(trt2), NA, val)) %>% mutate(trt1 = factor(df$trt1, levels = order$treat),
                                                                                              trt2 = factor(df$trt2, levels = rev(order$treat)),
                                                                                              sig = ifelse(is.na(val),1,ifelse(lwr < statsig & upr < statsig | lwr > statsig & upr > statsig,2,0)))  #Identify stat sig results
    ext <- data.frame(lwr=NA, med=NA, upr=NA, trt1=order$treat, trt2=order$treat, val=order$treat,sig = NA)
    
    df = rbind(df,ext) %>% slice(-1)
  }
  
  
  
  windows()
  ggplot(df, aes(x = trt1, y = trt2)) +
    geom_tile(aes(fill = sig), colour = ifelse(is.na(df$val),"white","black"), show.legend = F) +
    geom_text(aes(label = df$val), size = textsize) +
    ggtitle("League Table") +
    # add colour scale
    scale_fill_gradient2(name="",
                         # midpoint = 0.5,
                         na.value = colour,
                         low = "white",
                         mid = "white",
                         high = ifelse(!is.na(df$val),"lightgrey","white")) +
    # # move x-axis label to top
    # scale_x_discrete(position = "top") + xlab("NMA Model") +
    # use a white background, remove axis borders, labels, tickts
    theme_bw() + xlab(" ") + ylab(" ") +
    theme(panel.border=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  
}


###=========================NMA heatplot====

heatplot = function(data, border = "white", border_size = 2){
  
  #Arrange treatments according to SUCRA in primary analysis    
  data = data %>% arrange(-data[,2])  
  
  # Convert to long format, create labels for cells and factor variables in order of appearence (as_factor instead of as.factor). Aaron, I'm sure there is
  # a way to do this that doesn't require forcats? It's a Hadley package so I trust it to be maintained but base would be nicer
  
  long = data %>% gather(variable,value,-Treatment) %>% mutate(value2 = ifelse(!is.na(value),paste(round(value*100,0),"%",sep = ""),
                                                                               paste(round(value*100,0)))) %>% mutate(variable = as_factor(variable))
  
  #Orders treatments based on sucra. Rev needed because axis is flipped?
  long$Treatment = factor(data$Treatment,
                          levels = rev(data$Treatment))  
  
  
  
  long %>% ggplot(aes(y = Treatment, x = variable)) + 
    geom_tile(aes(fill = value), colour = "white", size = 2) + 
    geom_text(aes(label = value2)) + scale_fill_gradient2(name="Legend\n",
                                                          midpoint = 0.5,
                                                          limits = c(0, 1),
                                                          breaks = 0.5*0:2,
                                                          labels = percent(0.5*0:2),
                                                          na.value = I(rgb(255, 255, 255, maxColorValue=255)),
                                                          low = I(rgb(248, 105, 107, maxColorValue=255)),
                                                          mid = I(rgb(255, 235, 132, maxColorValue=255)),
                                                          high = I(rgb(0, 192, 82, maxColorValue=255))) +
    # move x-axis label to top
    scale_x_discrete(position = "top") + xlab("NMA Model") +
    # use a white background, remove borders
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.border=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    theme(legend.title =  element_text(face = "bold", size = 10),
          legend.text = element_text(size = 7.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
}

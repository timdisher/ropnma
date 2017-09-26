prep_wb = function(data, smd = FALSE,convert_contrast = TRUE,md = TRUE){
  data$trt_group = fct_infreq(data$trt_group) %>% droplevels(data)#reorders factor from most to least common
  
  data = droplevels(data) # Drops factor levels that don't exist (otherwise they are carried over)
  
  if(md == TRUE){
  if(convert_contrast == TRUE){
    if(smd == TRUE) long_wb_smd(data, hedge = TRUE) else long_wb(data = data)}else long_arm_wb(data = data)
  
  } else{
  
    long_wb(data = data, md = FALSE)
  }
}

nma = function(data,variable,drop,SA = FALSE, inc = FALSE, models = model){
  data = data %>% left_join(rop_data_study[c("studlab","design")],by = "studlab")
  
  if(SA == TRUE){
    data = data %>% filter(data[variable] != drop)} else data = data
    
    pw = pairwise(data = data,treat = trt_group, n = sample_size,mean = mean, sd = std_dev, studlab = studlab, sm = "MD")
    
    nc = netconnection(data = pw, treat1,treat2)
    
    if(nc$n.subnets > 1) return("This sensitivity analysis disconnects the network") else( data = prep_wb(data))
    
    data$xo = data$wide %>% left_join(rop_data_study %>% select(studlab,design), by = "studlab") %>%
      
      left_join(rop_data_arm %>% select(studlab,p_value) %>% distinct() %>% filter(!is.na(p_value)),by = "studlab") %>%
      
      mutate(se_2 = ifelse(design == "Crossover",se_paired(y_2,p_value,n_1),se_2)) %>%
      
      select(matches("t_"),matches("y_"),matches("se_"),V,na) %>% select(-y_1)
    
    
    if("y_4" %in% colnames(data$xo))  model = list(models$re,models$re_inc) else model = list(models$re3,models$re3_inc)
    
    if(inc == FALSE){
      nma = nma_cont(data, data$xo, data$treatments,params = params.re, model = model[[1]],bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 10, FE = FALSE, inc = inc)
      
    } else nma = nma_cont(data$wide, data$xo, data$treatments,params = params.re, model = model,bugsdir = bugsdir,n.iter = 100000, n.burnin = 40000,n.thin = 10, FE = FALSE, inc = inc)
    list(data = data,nma = nma)
}

#============================================================================================
#
#
#             ---- Arrange long data into winbugs format -----
#
#
#============================================================================================

#data = pa_reac, a data frame limited to the outcome of interest



long_wb = function(data = pa_reac,studlab = "studlab", trt = "trt_group",mean = "mean",sd = "std_dev",sample_size = "samplesize",
                   n_event = "num_events",md = TRUE){
  
  if(md == TRUE){
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
      
      winbugs_ready = wide_con %>% select(t_1:t_2,y_2,se_2,na)} 
    
  out = list(input = input,treatments = treatments,arm_wide = wide,wide = wide_con,wb = winbugs_ready)    
  } else{
       
    (input = data %>% select(studlab,trt_group,n_event,sample_size) %>% 
       rename(r = num_events,
              n = sample_size) %>%  
       arrange(studlab,trt_group) %>% ##Arranges treatments by factor level
       group_by(studlab) %>%
       mutate(arm = row_number(),
              t = as.numeric(trt_group),
              na = n()))
    
    treatments = input %>% ungroup() %>% select(trt_group,t) %>% distinct() %>% arrange(t) %>%
      rename(description = trt_group)
    
    (wide= input %>% select(-trt_group) %>%
        recast(studlab ~ variable + arm, id.var = c("studlab","arm")) %>% select(studlab:na_1) %>% rename(na = na_1)
    )
    
    wb = wide %>% select(-studlab)
    
    out = list(input = input,treatments = treatments,wide = wide, wb = wb) }
      
    
      
#####Output a list
  

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
  








#===================================================================================================









#=====================================================================
#
#
# Set up data for R2WinBUGS
#
#
#=====================================================================


#================
# Create data list
#===============

nma_winbugs_datalist = function(data,treatments,contrast = TRUE, md = TRUE){


#Find column indexes
if(md == TRUE){  
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

data = list(nt=nt,t=t, y=y, se=se, na=na, ns = ns)}
  } else{
    
    t_loc = grep("t_|t. != studlab",colnames(data),fixed = F)
    r_loc = grep("r_|r.",colnames(data),fixed = F)
    n_loc = grep("n_",colnames(data),fixed = F)
    na_loc = grep("na",colnames(data),fixed = F)
  
    t = as.matrix(data[,t_loc])
    
    # number of treatments
    nt = length(treatments$description)
    
    r = as.matrix(data[,r_loc])
    
    n = as.matrix(data[,n_loc])
    
    na = as.vector(data[,na_loc])
    
    ns = as.vector(length(r[,1]))
    
    data = list(nt = nt, ns = ns, r = r, n = n, na = na,t = t)
      
  }
    
    
  data
  
}





#============================================================================================
#
#
#                                 Run the analysis 
#
# - Inconsistency analysis not coded for datasets with no four arm trials
#============================================================================================


nma_cont = function(names,data,treatments,n.iter = 40000, n.burnin = 20000, model, params, 
                    FE = TRUE,inc = FALSE, bugsdir = "c:/Users/TheTimbot/Desktop/WinBUGS14", n.thin = 1, debug = F, cont = TRUE){

  if(cont == TRUE){
  data = nma_winbugs_datalist(data,treatments)} else data = nma_winbugs_datalist(data,treatments, md = FALSE)
  
  model_con = model[[1]]
  
  
  model.con = bugs(data, NULL, params, model.file= model_con,
               n.chains = 3, n.iter = n.iter, n.burnin = n.burnin, n.thin= n.thin, 
               bugs.directory = bugsdir, debug= debug)
  
  if(inc == TRUE){
    model_inc = model[[2]]
    params_inc = c("sd","resdev","dev","totresdev","d")
  model.inc = bugs(data, NULL, params_inc, model.file= model_inc,
                   n.chains = 3, n.iter = n.iter, n.burnin = n.burnin, n.thin= n.thin, 
                   bugs.directory = bugsdir, debug=F) }



#============================================================
# Create outputs
#============================================================

#======================================
# Output all pairwise comparisons
#======================================
  
  if(cont == TRUE){
# summarize mean differences - mean, 2.5%, median, 97.5%
  md = (as.matrix(model.con$summary[grep("meandif", rownames(model.con$summary), fixed=F), c(1,3,5,7)]))} else{
    
    md = (as.matrix(model.con$summary[grep("or", rownames(model.con$summary), fixed=F), c(1,3,5,7)])) 
  }
  
# summarize the probability better - mean, 2.5%, median, 97.5%
  better = (as.matrix(model.con$summary[grep("better", rownames(model.con$summary), fixed=F), c(1,2)]))
  
  
# conver comparison columns into a matrix
  if(cont == TRUE){
  b <- gsub("\\meandif|\\[|\\]", "",rownames(md))} else b <- gsub("\\or|\\[|\\]", "",rownames(md))
  b_clean <- matrix(as.numeric(unlist(strsplit(b, ","))), ncol=2, byrow=T)
  comp <- paste(treatments$description[b_clean[,2]], ' vs.', treatments$description[b_clean[,1]])
  
  md_med <- round(md[,"50%"],2)
  
  md_ci <- paste(paste(round(md[,'2.5%'],2), round(md[,'97.5%'],2), sep =" to "), sep="")
  
  hr_med_ci = paste(round(md[,"50%"],2) , paste("(", paste(round(md[,'2.5%'],2), round(md[,'97.5%'],2), sep =" to "), ")", sep=""))
  
  better_mean_sd = paste(round(better[,"mean"],2) , paste("(", paste(round(better[,'sd'],2), ")", sep=""), sep=""))
  
  better_mean = round(better[,"mean"],2)
  
  better_sd = round(better[,'sd'],2)
  
  if(cont == TRUE){
  comp = cbind(comp, md_med, md_ci, better_mean, better_sd); 
  colnames(comp) = c('Comparison (Trt A vs. Trt B)', 'Mean Difference of Trt A vs. Trt B', "95% CrI of Mean Difference", 'Probability of Trt A better than Trt B', "SD of Probability Better")} else{
    
    comp = cbind(comp, md_med, md_ci, better_mean, better_sd); colnames(comp) = c('Comparison (Trt A vs. Trt B)', 'OR of Trt A vs. Trt B', "95% CrI of OR", 'Probability of Trt A better than Trt B', "SD of Probability Better")  
  }
  rownames(comp) <- NULL


#======================================
#SUCRA, Rankings, and Probability Best Plots
#======================================

  
# sucra summary
sucra = model.con$summary[grep("SUCRA", rownames(model.con$summary), fixed=F),] #mean SUCRA, remove ",1" to get median values

# ranks summary
ranks = model.con$summary[grep("rk", rownames(model.con$summary), fixed=F),] #mean rank, remove ",1" to get median values

# prob best summary
probs = model.con$summary[grep("best", rownames(model.con$summary), fixed=F),] #mean probability best, remove ",1" to get median values

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
probs <- model.con$mean$prob; rownames(probs) <- treatments$description; colnames(probs) <- seq(1,data$nt)

dat <- melt(cbind(probs)); colnames(dat) <- c('Treatments', 'Rankings', 'Probability of Best Treatment')


rankogram = ggplot(dat,aes(x = Rankings, y = `Probability of Best Treatment`,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format()) + scale_x_continuous(breaks=seq(1,data$nt))




#=====================================================
# Inconsistency results
#=====================================================
if(inc == TRUE){
  
  if(cont == TRUE){
  inc <- model.inc$summary[grep("d[", rownames(model.inc$summary), fixed=T), c(1,3,5,7)]} else inc <- exp(model.inc$summary[grep("d[", rownames(model.inc$summary), fixed=T), c(1,3,5,7)]) 

  inc <- cbind(inc, md,comp)

  d<-inc[!(as.numeric(inc[,4])>60),]

  inc_ci = paste(round(as.numeric(d[,3]),2), paste("(", paste(round(as.numeric(d[,2]),2), round(as.numeric(d[,4]),2), sep =" to "), ")", sep=""))
  con_ci = paste(round(as.numeric(d[,7]),2), paste("(", paste(round(as.numeric(d[,6]),2), round(as.numeric(d[,8]),2), sep =" to "), ")", sep=""))

  inc.con = cbind(con_ci, inc_ci)
  inc.con = rbind(inc.con, c(paste('DIC=', round(model.con$DIC,2), ' totresdev=', round(model.con$summary['totresdev',1],2)), 
                             paste('DIC=', round(model.inc$DIC,2), ' totresdev=', round(model.inc$summary['totresdev',1],2))))
  rownames(inc.con) = c(d[,9], 'Model Fit Statistics')
  colnames(inc.con) = c('FE Consistency NMA', 'FE Inconsistency NMA')




  x = model.con$summary[grep("^dev[[]", rownames(model.con$summary), fixed=F),1];

  
  
  y = model.inc$summary[grep("^dev[[]", rownames(model.inc$summary), fixed=F),1]
  
  names$V1 = c(1:length(names$studlab))
  
  treats = names %>% select(studlab,matches("t_")) %>% select(-t_1)   
  
  
  
  if("t_4" %in% colnames(treats)){
  treats["t_2"] = treatments[match(treats[["t_2"]], treatments[["t"]]),1]
  treats["t_3"] = treatments[match(treats[["t_3"]], treatments[["t"]]),1]
  treats["t_4"] = treatments[match(treats[["t_4"]], treatments[["t"]]),1]
  
  dict = treats %>% rename('2' = t_2,
                           '3' = t_3,
                           '4' = t_4) %>%gather(arm,treat,2:4) %>% mutate(arm = as.numeric(arm)) %>% rename(V2 = arm)} else{
                             
                             treats["t_2"] = treatments[match(treats[["t_2"]], treatments[["t"]]),1]
                             treats["t_3"] = treatments[match(treats[["t_3"]], treatments[["t"]]),1]
                             
                             dict = treats %>% rename('2' = t_2,
                                                      '3' = t_3) %>%gather(arm,treat,2:3) %>% mutate(arm = as.numeric(arm)) %>% rename(V2 = arm)                      
                             
                             
                           }
  

  ric = as.data.frame(cbind(x,y)) %>% rownames_to_column("Arm") %>% rename('Consistency NMA' = x,
                                                                           'Inconsistency NMA' = y) 
  ric = ric %>%
 bind_cols(as.data.frame(str_extract_all(ric$Arm, regex('(?<=\\[|,)[0-9]+(?=\\]|-?)'),simplify = TRUE))) %>% 
    mutate(V1 = as.integer(as.character(V1)),V2 = as.integer(as.character(V2))) %>%
    left_join(names %>% select(studlab, V1), by = c("V1")) %>% left_join(dict, by = c("studlab","V2")) %>%
    unite(label,studlab,treat,sep = ", ")





  p <- ggplot(as.data.frame(ric), aes(`Consistency NMA`, `Inconsistency NMA`), label= label)



  #Scatterplot for RE Model

  dev_plot = p + geom_point() + scale_y_continuous(limits=c(0, 5)) + scale_x_continuous(limits=c(0, 5)) +  geom_abline(intercept = 0, colour='gray56', lty=2) +
    ggtitle("Deviance Plot") +  theme(plot.title = element_text(hjust = 0.5)) +
    geom_text(aes(label=ifelse(`Consistency NMA`>1.25|`Consistency NMA`<0.25,as.character(label),'')),hjust=0, nudge_x = 0.05, angle=0, size=2.5)
  
  
  if(FE == TRUE){  
    results = list(model = model,comp = comp,rr = rr,rankogram = rankogram)
    
  }else {sd = model.con$summary[grep("^sd$",row.names(model.con$summary),fixed = F),]
  
  results = list(model = model,sd = sd,comp = comp,rr = rr,rankogram = rankogram, inc_table = inc.con,
                 devplot = dev_plot)}
  
  
  results
} else if(FE == TRUE){  
    results = list(model = model,comp = comp,rr = rr,rankogram = rankogram)
    
  }else {sd = model.con$summary[grep("^sd$",row.names(model.con$summary),fixed = F),]
  
  results = list(model = model,sd = sd,comp = comp,rr = rr,rankogram = rankogram)}


results
}

  

#=====================================================================================


#                                     -Models -


#=====================================================================================

fe_binom = function(){
  
  # *** PROGRAM STARTS
  for(i in 1:ns){                 # LOOP THROUGH STUDIES
    mu[i] ~ dnorm(0,.0001)      # vague priors for all trial baselines
    for (k in 1:na[i])  {       # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])    # binomial likelihood
      # model for linear predictor
      logit(p[i,k]) <- mu[i] + d[t[i,k]] - d[t[i,1]]
      # expected value of the numerators 
      rhat[i,k] <- p[i,k] * n[i,k]
      #Deviance contribution
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))
                       +  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))
    }
    # summed residual deviance contribution for this trial
    resdev[i] <- sum(dev[i,1:na[i]])
  }   
  totresdev <- sum(resdev[])      # Total Residual Deviance
  d[1]<-0    # treatment effect is zero for reference treatment
  # vague priors for treatment effects
  for(k in 2:nt){  d[k] ~ dnorm(0,.0001) }
  
  ################################################################################ 
  # Extra code for all odds ratios and log odds ratios, ranking, and absolute effects, and relative effects 
  # on alternative scales: Numbers Needed to Treat, Risk Difference, Relative Risks 
  ################################################################################ 
  # pairwise ORs and LORs for all possible pair-wise comparisons, if nt>2 
  for(c in 1:(nt-1)) { for (k in (c+1):nt)
  {
    or[c,k] <- exp(d[k] - d[c]) 
    lor[c,k] <- (d[k]-d[c])
    better[c,k] <- 1 - step(d[k] - d[c])
  } }

  # ranking on relative scale 
  for (k in 1:nt) { 
  #rk[k] <- nt+1-rank(d[],k) # assumes events are “good”
  rk[k] <- rank(d[],k) # assumes events are “bad”
  
  best[k] <- equals(rk[k],1) #calculate probability that treat k is best 
  for (h in 1:nt){ prob[h,k] <- equals(rk[k],h)
  } # calculates probability that treat k is h-th best 
  }
  
  for (k in 1:nt) {
    for (h in 1:nt) {
      cumeffectiveness[k, h] <- sum(prob[k, 1:h])
    }
  }
  for (i in 1:nt) {
    SUCRA[i] <- sum(cumeffectiveness[i, 1:(nt - 1)])/(nt - 
                                                        1)
  }


}


re_binom = function(){
  # *** PROGRAM STARTS
  for(i in 1:ns){                      # LOOP THROUGH STUDIES
    w[i,1] <- 0    # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)           # vague priors for all trial baselines
    for (k in 1:na[i]) {             # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]  # model for linear predictor
      rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
      #Deviance contribution
      dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))  
                       +  (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))         }
    #  summed residual deviance contribution for this trial
    resdev[i] <- sum(dev[i,1:na[i]])       
    for (k in 2:na[i]) {             # LOOP THROUGH ARMS
      # trial-specific LOR distributions
      delta[i,k] ~ dnorm(md[i,k],taud[i,k])
      # mean of LOR distributions (with multi-arm trial correction)
      md[i,k] <-  d[t[i,k]] - d[t[i,1]] + sw[i,k]
      # precision of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k
      # adjustment for multi-arm RCTs
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])
      # cumulative adjustment for multi-arm trials
      sw[i,k] <- sum(w[i,1:k-1])/(k-1)
    }
  }   
  totresdev <- sum(resdev[])           # Total Residual Deviance
  d[1]<-0       # treatment effect is zero for reference treatment
  # vague priors for treatment effects
  for (k in 2:nt){  d[k] ~ dnorm(0,.0001) }
  sd ~ dunif(0,5)     # vague prior for between-trial SD
  tau <- pow(sd,-2)   # between-trial precision = (1/between-trial variance)
  
  ################################################################################ 
  # Extra code for all odds ratios and log odds ratios, ranking, and absolute effects, and relative effects 
  # on alternative scales: Numbers Needed to Treat, Risk Difference, Relative Risks 
  ################################################################################ 
  # pairwise ORs and LORs for all possible pair-wise comparisons, if nt>2 
  for(c in 1:(nt-1)) { for (k in (c+1):nt)
  {
    or[c,k] <- exp(d[k] - d[c]) 
    lor[c,k] <- (d[k]-d[c])
    better[c,k] <- 1 - step(d[k] - d[c])
  } }
  
  # ranking on relative scale 
  for (k in 1:nt) { 
    #rk[k] <- nt+1-rank(d[],k) # assumes events are “good”
    rk[k] <- rank(d[],k) # assumes events are “bad”
    
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best 
    for (h in 1:nt){ prob[h,k] <- equals(rk[k],h)
    } # calculates probability that treat k is h-th best 
  }
  
  for (k in 1:nt) {
    for (h in 1:nt) {
      cumeffectiveness[k, h] <- sum(prob[k, 1:h])
    }
  }
  for (i in 1:nt) {
    SUCRA[i] <- sum(cumeffectiveness[i, 1:(nt - 1)])/(nt - 
                                                        1)
  }

}


#==============================
# Write model files
#==============================
normal_models = function(fe = fe_binom, re = re_binom){
  
  write.model(fe, "fe-binom.txt")
  MODELFILE.fe <- c("fe-binom.txt")
  
  write.model(re, "re-binom.txt")
  MODELFILE.re <- c("re-binom.txt")
  

  
  list = list(fe = MODELFILE.fe,re = MODELFILE.re )
  
  list
}

#==============================
#=============================================================================================


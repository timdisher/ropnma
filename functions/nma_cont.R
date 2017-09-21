#Functions for NMA


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
#
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

#============================================================================================
#
#
#             ---- Arrange long data into winbugs format -----
#
#
#============================================================================================

#data = pa_reac, a data frame limited to the outcome of interest



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

long_wb_smd= function(data = pa_reac,studlab = "studlab", trt = "trt_group",mean = "mean",sd = "std_dev",sample_size = "samplesize"){
  
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
                                                        
                                                        (wide = wide %>%  mutate(se_2 = se_smd(y_2,sd_2,sd_1,n_2,n_1) %>% arrange(na)))}

    
     #arrange by number of arms (see WiBUGS code for why)
  
  wide %>% select(studlab,matches("t_"),matches("y_"),matches("se_"),V,na) %>% select(-y_1)
  

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

  if(ns4 >0)  data = list(nt=nt, ns2=ns2, ns3=ns3, ns4=ns4,t=t, y=y, se=se, na=na, V=V) else data = list(nt=nt, ns2=ns2, ns3=ns3,t=t, y=y, se=se, na=na, V=V)
  
 
  
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





#============================================================================================
#
#
#                                 Run the analysis 
#
# - Inconsistency analysis not coded for datasets with no four arm trials
#============================================================================================


nma_cont = function(names,data,treatments,n.iter = 40000, n.burnin = 20000, model, params, 
                    FE = TRUE,inc = FALSE, bugsdir = "c:/Users/TheTimbot/Desktop/WinBUGS14", n.thin = 1){

  data = nma_winbugs_datalist(data,treatments)
  
  model_con = model[[1]]
  
  
  model.con = bugs(data, NULL, params, model.file= model_con,
               n.chains = 3, n.iter = n.iter, n.burnin = n.burnin, n.thin= n.thin, 
               bugs.directory = bugsdir, debug=F)
  
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
  
# summarize mean differences - mean, 2.5%, median, 97.5%
  md = (as.matrix(model.con$summary[grep("meandif", rownames(model.con$summary), fixed=F), c(1,3,5,7)]))
  
# summarize the probability better - mean, 2.5%, median, 97.5%
  better = (as.matrix(model.con$summary[grep("better", rownames(model.con$summary), fixed=F), c(1,2)]))
  
  
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
  
  inc <- model.inc$summary[grep("d[", rownames(model.inc$summary), fixed=T), c(1,3,5,7)]

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


#                                     -Models - 3 arms-


#=====================================================================================

#=====================================================
# Normal likelihood, gaussian link 
# Fixed effects model for multi-arm trials
#=====================================================


fe_normal_gaus_3arm <- function()	                             
  # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol
  # Website: www.bris.ac.uk/cobm/research/mpes					
{								
  # *** PROGRAM STARTS								
  # LOOP THROUGH 2-ARM STUDIES								
  for(i in 1:ns2) { 								
    
    # normal likelihood for 2-arm trials	
    y[i,2] ~ dnorm(delta[i,2],prec[i,2])
    # Deviance contribution for trial i			
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]      #TSD2 - pg 20, table				
  }
  
  #LOOP THROUGH 3-ARM STUDIES								
  for(i in (ns2+1):(ns2+ns3)) {								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,(k+1)]*equals(j,k) }       # TSD 2 - pg 37, distribution of Y_i,xxx 				
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])	#TSD2 - pg 20, table							
  }
  
  
  #LOOP THROUGH ALL STUDIES								
  for(i in 1:(ns2+ns3)){								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      delta[i,k] <- d[t[i,k]] - d[t[i,1]]								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$							
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$							
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k] <- 1 - step(d[k] - d[c])
      
    }								
  }								
  
  # ranking calculations								
  for (k in 1:nt){								
    rk[k] <- rank(d[],k) # assumes differences < 0 favor the comparator				
    
    # Prob Best
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best  			
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h)								
    }								
  }	
  
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program							


#=====================================================
# Normal likelihood, gaussian link 
# Fixed effects model for multi-arm trials inconsistency
#=====================================================
fe_normal_gaus__3arminc <- function(){
  for(i in 1:ns2) {                    # LOOP THROUGH 2-ARM STUDIES
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials
    #Deviance contribution for trial i
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]
  }
  for(i in (ns2+1):(ns2+ns3)) {        # LOOP THROUGH THREE-ARM STUDIES
    for (k in 1:(na[i]-1)) {    # set variance-covariance matrix
      for (j in 1:(na[i]-1)) {
        Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k)
      }
    }
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])  #Precision matrix
    # multivariate normal likelihood for 3-arm trials   
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)]) 
    #Deviance contribution for trial i
    for (k in 1:(na[i]-1)){  # multiply vector & matrix
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])
    }
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])
  }
  
 
  for(i in 1:(ns2+ns3)){                      #   LOOP THROUGH ALL STUDIES
    for (k in 2:na[i]) {             #  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2)   # calculate variances
      prec[i,k] <- 1/var[i,k]      # set precisions
      # trial-specific LOR distributions
      delta[i,k] <- d[t[i,1],t[i,k]]
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]
    }
  }   
  totresdev <- sum(resdev[])            #Total Residual Deviance
  for (k in 1:nt) { d[k,k] <- 0 }
  for (c in 1:(nt-1)) {  # priors for all mean treatment effects
    for (k in (c+1):nt)  { 
      d[c,k] ~ dnorm(0,.001) 
      hr[c,k] <- exp(d[c,k])
    } 
  }  
  
}                                     # *** PROGRAM ENDS 




#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials
#=====================================================


re_normal_gaus_3arm <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  
  for(i in 1:(ns2+ns3)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects								
  sd ~ dunif(0,5) # vague prior for between-trial SD								
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program							



#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials using informative priors on SMD
#=====================================================


re_normal_gaus_3arm_inform <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  for(i in 1:(ns2+ns3)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  prior.prec <- 1/(2.41*2.41)
  sd ~ dt(-4.52, prior.prec,5) # vague prior for between-trial SD	
  tausq <- exp(sd)
  tau <- pow(tausq,-1)# between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program	


#====================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials inconsistency
#====================================
re_normal_gaus_3arm_inc <- function(){
  for(i in 1:ns2) {                    # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    #Deviance contribution for trial i								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]								
  }								
  for(i in (ns2+1):(ns2+ns3)) {        # LOOP THROUGH THREE-ARM STUDIES								
    for (k in 1:(na[i]-1)) {    # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) {								
        Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k)								
      }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])  #Precision matrix								
    # multivariate normal likelihood for 3-arm trials   								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)]) 								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){  # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }						
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  
  for(i in 1:(ns2+ns3)){                      #   LOOP THROUGH ALL STUDIES								
    for (k in 2:na[i]) {             #  LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2)   # calculate variances								
      prec[i,k] <- 1/var[i,k]      # set precisions								
      # trial-specific LOR distributions								
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
  }   								
  totresdev <- sum(resdev[])            #Total Residual Deviance								
  for (k in 1:nt) { d[k,k] <- 0 }								
  for (c in 1:(nt-1)) {  # priors for all mean treatment effects								
    for (k in (c+1):nt)  { 
      d[c,k] ~ dnorm(0,.001) 
      hr[c,k] <- exp(d[c,k])
    } 								
  }  								
  sd ~ dunif(0,5)     # vague prior for between-trial SD								
  tau <- pow(sd,-2)   # between-trial precision = (1/between-trial variance)								
}                                     # *** PROGRAM ENDS




#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials with continuous meta-regression
#=====================================================


re_normal_gaus_metareg_3arm <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix	
    
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  
  for(i in 1:(ns2+ns3)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001)
    beta[k] <- B} # vague priors for treatment effects
  B ~ dnorm(0,.0001)
  sd ~ dunif(0,5) # vague prior for between-trial SD								
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program	


#=====================================================================================


#                                     -Models - 4 arms-


#=====================================================================================

#=====================================================
# Normal likelihood, gaussian link 
# Fixed effects model for multi-arm trials
#=====================================================


fe_normal_gaus <- function()	                             
  # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol
  # Website: www.bris.ac.uk/cobm/research/mpes					
{								
  # *** PROGRAM STARTS								
  # LOOP THROUGH 2-ARM STUDIES								
  for(i in 1:ns2) { 								
    
    # normal likelihood for 2-arm trials	
    y[i,2] ~ dnorm(delta[i,2],prec[i,2])
    # Deviance contribution for trial i			
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]      #TSD2 - pg 20, table				
  }
  
  #LOOP THROUGH 3-ARM STUDIES								
  for(i in (ns2+1):(ns2+ns3)) {								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,(k+1)]*equals(j,k) }       # TSD 2 - pg 37, distribution of Y_i,xxx 				
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])	#TSD2 - pg 20, table							
  }
  
  #LOOP THROUGH 4-ARM STUDIES								
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) {								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,(k+1)]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }	
  
  #LOOP THROUGH ALL STUDIES								
  for(i in 1:(ns2+ns3+ns4)){								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      delta[i,k] <- d[t[i,k]] - d[t[i,1]]								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$							
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$							
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k] <- 1 - step(d[k] - d[c])
      
    }								
  }								
  
  # ranking calculations								
  for (k in 1:nt){								
    rk[k] <- rank(d[],k) # assumes differences < 0 favor the comparator				
    
    # Prob Best
    best[k] <- equals(rk[k],1) #calculate probability that treat k is best  			
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h)								
    }								
  }	
  
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program							


#=====================================================
# Normal likelihood, gaussian link 
# Fixed effects model for multi-arm trials inconsistency
#=====================================================
fe_normal_gaus_inc <- function(){
  for(i in 1:ns2) {                    # LOOP THROUGH 2-ARM STUDIES
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials
    #Deviance contribution for trial i
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]
  }
  for(i in (ns2+1):(ns2+ns3)) {        # LOOP THROUGH THREE-ARM STUDIES
    for (k in 1:(na[i]-1)) {    # set variance-covariance matrix
      for (j in 1:(na[i]-1)) {
        Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k)
      }
    }
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])  #Precision matrix
    # multivariate normal likelihood for 3-arm trials   
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)]) 
    #Deviance contribution for trial i
    for (k in 1:(na[i]-1)){  # multiply vector & matrix
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])
    }
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])
  }
  
  #LOOP THROUGH 4-ARM STUDIES								
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) {								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }	
  
  for(i in 1:(ns2+ns3+ns4)){                      #   LOOP THROUGH ALL STUDIES
    for (k in 2:na[i]) {             #  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2)   # calculate variances
      prec[i,k] <- 1/var[i,k]      # set precisions
      # trial-specific LOR distributions
      delta[i,k] <- d[t[i,1],t[i,k]]
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]
    }
  }   
  totresdev <- sum(resdev[])            #Total Residual Deviance
  for (k in 1:nt) { d[k,k] <- 0 }
  for (c in 1:(nt-1)) {  # priors for all mean treatment effects
    for (k in (c+1):nt)  { 
      d[c,k] ~ dnorm(0,.001) 
      hr[c,k] <- exp(d[c,k])
    } 
  }  
  
}                                     # *** PROGRAM ENDS 




#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials
#=====================================================


re_normal_gaus <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) { # LOOP THROUGH 4-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  
  
  for(i in 1:(ns2+ns3+ns4)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects								
  sd ~ dunif(0,5) # vague prior for between-trial SD								
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program							



#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials using informative priors on SMD
#=====================================================


re_normal_gaus_inform <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix								
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) { # LOOP THROUGH 4-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  
  
  for(i in 1:(ns2+ns3+ns4)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment								
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  prior.prec <- 1/(2.41*2.41)
  sd ~ dt(-4.52, prior.prec,5) # vague prior for between-trial SD	
  tausq <- exp(sd)
  tau <- pow(tausq,-1)# between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program	


#====================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials inconsistency
#====================================
re_normal_gaus_inc <- function(){
  for(i in 1:ns2) {                    # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    #Deviance contribution for trial i								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2]								
  }								
  for(i in (ns2+1):(ns2+ns3)) {        # LOOP THROUGH THREE-ARM STUDIES								
    for (k in 1:(na[i]-1)) {    # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) {								
        Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k)								
      }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])  #Precision matrix								
    # multivariate normal likelihood for 3-arm trials   								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)]) 								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){  # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }						
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) { # LOOP THROUGH 4-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  for(i in 1:(ns2+ns3+ns4)){                      #   LOOP THROUGH ALL STUDIES								
    for (k in 2:na[i]) {             #  LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2)   # calculate variances								
      prec[i,k] <- 1/var[i,k]      # set precisions								
      # trial-specific LOR distributions								
      delta[i,k] ~ dnorm(d[t[i,1],t[i,k]] ,tau)
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
  }   								
  totresdev <- sum(resdev[])            #Total Residual Deviance								
  for (k in 1:nt) { d[k,k] <- 0 }								
  for (c in 1:(nt-1)) {  # priors for all mean treatment effects								
    for (k in (c+1):nt)  { 
      d[c,k] ~ dnorm(0,.001) 
   
    } 								
  }  								
  sd ~ dunif(0,5)     # vague prior for between-trial SD								
  tau <- pow(sd,-2)   # between-trial precision = (1/between-trial variance)								
}                                     # *** PROGRAM ENDS

#=====================================================
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials with continuous meta-regression
#=====================================================


re_normal_gaus_metareg <- function()	                             # this code for this model was adapted from WinBUGS code from the multi-parameter Evidence Synthesis Research Group at the University of Bristol:  Website: www.bris.ac.uk/cobm/research/mpes					
{								
  for(i in 1:ns2) { # LOOP THROUGH 2-ARM STUDIES								
    y[i,2] ~ dnorm(delta[i,2],prec[i,2]) # normal likelihood for 2-arm trials								
    resdev[i] <- (y[i,2]-delta[i,2])*(y[i,2]-delta[i,2])*prec[i,2] #Deviance contribution for trial i								
  }								
  for(i in (ns2+1):(ns2+ns3)) { # LOOP THROUGH 3-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,]) #Precision matrix	
    
    # multivariate normal likelihood for 3-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega[i,1:(na[i]-1),1:(na[i]-1)])
    
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }		
  
  for(i in (ns2+ns3+1):(ns2+ns3+ns4)) { # LOOP THROUGH 4-ARM STUDIES								
    for (k in 1:(na[i]-1)) { # set variance-covariance matrix								
      for (j in 1:(na[i]-1)) { Sigma2[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }								
    }								
    Omega2[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma2[i,,]) #Precision matrix								
    # multivariate normal likelihood for 4-arm trials								
    y[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]],Omega2[i,1:(na[i]-1),1:(na[i]-1)])								
    #Deviance contribution for trial i								
    for (k in 1:(na[i]-1)){ # multiply vector & matrix								
      ydiff[i,k]<- y[i,(k+1)] - delta[i,(k+1)]								
      z[i,k]<- inprod2(Omega2[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])								
    }								
    resdev[i]<- inprod2(ydiff[i,1:(na[i]-1)], z[i,1:(na[i]-1)])								
  }
  
  
  
  for(i in 1:(ns2+ns3+ns4)){ # LOOP THROUGH ALL STUDIES								
    w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm								
    delta[i,1] <- 0 # treatment effect is zero for control arm								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      var[i,k] <- pow(se[i,k],2) # calculate variances								
      prec[i,k] <- 1/var[i,k] # set precisions								
      dev[i,k] <- (y[i,k]-delta[i,k])*(y[i,k]-delta[i,k])*prec[i,k]								
    }								
    #resdev[i] <- sum(dev[i,2:na[i]]) # summed residual deviance contribution for this trial								
    for (k in 2:na[i]) { # LOOP THROUGH ARMS								
      delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific treat effects distributions								
      md[i,k] <- d[t[i,k]] - d[t[i,1]] + (beta[t[i,k]]-beta[t[i,1]])*(x[i]-mx) + sw[i,k] # mean of treat effects distributions (with multi-arm trial correction)								
      taud[i,k] <- tau *2*(k-1)/k # precision of treat effects distributions (with multi-arm trial correction)								
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs								
      sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials								
    }								
  }								
  totresdev <- sum(resdev[]) #Total Residual Deviance								
  d[1]<-0 # treatment effect is zero for reference treatment
  beta[1] <- 0
  for (k in 2:nt){ d[k] ~ dnorm(0,.0001)
    beta[k] <- B} # vague priors for treatment effects
  B ~ dnorm(0,.0001)
  sd ~ dunif(0,5) # vague prior for between-trial SD								
  tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance)								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}	                                                        # END Program	



#=====================================================
# Normal likelihood, gaussian link- Data by arm
# Random effects model for multi-arm trials with baseline risk meta-regression
#=====================================================


re_normal_armdata_meta = function(){
  for(i in 1:ns){                      #   LOOP THROUGH STUDIES
    w[i,1] <- 0    # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.001)           # vague priors for all trial baselines
    for (k in 1:na[i]) {             #  LOOP THROUGH ARMS
      var[i,k] <- pow(se[i,k],2)   # calculate variances
      prec[i,k] <- 1/var[i,k]      # set precisions
      y[i,k] ~ dnorm(theta[i,k],prec[i,k]) # binomial likelihood
      theta[i,k] <- mu[i] + delta[i,k] + (beta[t[i,k]]-beta[t[i,1]])*(mu[i]-mx)  # model for linear predictor
      #Deviance contribution
      dev[i,k] <- (y[i,k]-theta[i,k])*(y[i,k]-theta[i,k])*prec[i,k]
    }
    #  summed residual deviance contribution for this trial
    resdev[i] <- sum(dev[i,1:na[i]])       
    for (k in 2:na[i]) {             # LOOP THROUGH ARMS
      # trial-specific LOR distributions
      delta[i,k] ~ dnorm(md[i,k],taud[i,k])
      # mean of LOR distributions, with multi-arm trial correction
      md[i,k] <-  d[t[i,k]] - d[t[i,1]] + sw[i,k]
      # precision of LOR distributions (with multi-arm trial correction)
      taud[i,k] <- tau *2*(k-1)/k
      # adjustment, multi-arm RCTs
      w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])
      # cumulative adjustment for multi-arm trials
      sw[i,k] <- sum(w[i,1:k-1])/(k-1)
    }
  }   
  totresdev <- sum(resdev[])            #Total Residual Deviance
  d[1]<-0       # treatment effect is zero for control arm
  beta[1] <- 0
  # vague priors for treatment effects
  for (k in 2:nt){  d[k] ~ dnorm(0,.0001)
    beta[k] <- B}
  sd ~ dunif(0,5)     # vague prior for between-trial SD
  B ~ dnorm(0,.0001)
  tau <- pow(sd,-2)   # between-trial precision = (1/between-trial variance)
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  # Extra code for all mean differences, rankings								
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$								
  
  # pairwise mean differences for all possible pair-wise comparisons, if nt>2								
  for (c in 1:(nt-1)) {								
    for (k in (c+1):nt) {								
      meandif[c,k] <- (d[k] - d[c])
      # pairwise comparison between all treatments
      better[c,k]  <- 1 - step(d[k] - d[c])
      
    }								
  }								
  # ranking calculations								
  for (k in 1:nt) {
    # assumes differences<0 favor the comparator   ===> number of elements in d[] that are less than or equal to d[k] 					
    rk[k] <- rank(d[],k) 			
    
    #calculate probability that treat k is best		===> k is the best treatment if rk[k] = 1
    best[k] <- equals(rk[k],1)
    
    for(h in 1:nt) {								
      prob[k,h]<- equals(rk[k],h) 				
    }								
  }								
  for(k in 1:nt) {								
    for(h in 1:nt) {								
      cumeffectiveness[k,h]<- sum(prob[k,1:h])								
    }								
  }								
  #SUCRAS#								
  for(i in 1:nt) {								
    SUCRA[i]<- sum(cumeffectiveness[i,1:(nt-1)]) /(nt-1)								
  }
}


#==============================
# Write model files
#==============================
normal_models = function(fe = fe_normal_gaus, re = re_normal_gaus, re_inf = re_normal_gaus_inform, 
                         fe_inc = fe_normal_gaus_inc, re_inc = re_normal_gaus_inc, re3 = re_normal_gaus_3arm, re3_inc = re_normal_gaus_3arm_inc,
                         re_meta = re_normal_gaus_metareg, re3_meta = re_normal_gaus_metareg_3arm, re_arm_meta = re_normal_armdata_meta){
  
  write.model(fe, "fe-normal-gaus.txt")
  MODELFILE.fe <- c("fe-normal-gaus.txt")
  
  write.model(re, "re-normal-gaus.txt")
  MODELFILE.re <- c("re-normal-gaus.txt")
  
  write.model(fe_inc, "fe-normal-gaus-inc.txt")
  MODELFILE.fe_inc <- c("fe-normal-gaus-inc.txt")
  
  
  write.model(re_inc, "re-normal-gaus-inc.txt")
  MODELFILE.re_inc <- c("re-normal-gaus-inc.txt")
  
  
  write.model(re_inf, "re-normal-gaus-inf.txt")
  MODELFILE.re_inf <- c("re-normal-gaus-inf.txt")
  
  write.model(re3, "re_normal_gaus_3arm.txt")
  MODELFILE.re3 <- c("re_normal_gaus_3arm.txt")
  
  
  write.model(re3_inc, "re_normal_gaus_3arm_inc.txt")
  MODELFILE.re3_inc <- c("re_normal_gaus_3arm_inc.txt")
  
  write.model(re_meta, "re_normal_metareg.txt")
  MODELFILE.re_meta <- c("re_normal_metareg.txt")
  
  write.model(re3_meta, "re_normal_metareg_3arm.txt")
  MODELFILE.re3_meta <- c("re_normal_metareg_3arm.txt")
  
  write.model(re_arm_meta, "re_normal_armdata_meta.txt")
  MODELFILE.re_armdata_meta <- c("re_normal_armdata_meta.txt")
  
  
  list = list(fe = MODELFILE.fe, fe_inc= MODELFILE.fe_inc,re = MODELFILE.re, 
              re_inf = MODELFILE.re_inf ,re_inc = MODELFILE.re_inc,re3 = MODELFILE.re3,
              re3_inc = MODELFILE.re3_inc, re_meta = MODELFILE.re_meta, re3_meta = MODELFILE.re3_meta, re_arm_meta =  MODELFILE.re_armdata_meta )
  
  list
}

#==============================
#=============================================================================================




#========================================
# Test zone
#=========================================
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
  results = list(model = model,comp = comp,rr = rr,rankogram = rankogram)
  
}else {sd = model$summary[grep("^sd$",row.names(model$summary),fixed = F),]

results = list(model = model,sd = sd,comp = comp,rr = rr,rankogram = rankogram)}


results
}

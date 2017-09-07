#Functions for NMA
data = pa_reac_wb


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#                              -----Models-----


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Normal likelihood, gaussian link 
# Fixed effects model for multi-arm trials
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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
      for (j in 1:(na[i]-1)) { Sigma[i,j,k] <- V[i]*(1-equals(j,k)) + var[i,k+1]*equals(j,k) }       # TSD 2 - pg 37, distribution of Y_i,xxx 				
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

write.model(fe_normal_gaus, "fe-normal-gaus.txt")
MODELFILE.fe <- c("fe-normal-gaus.txt")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Normal likelihood, gaussian link 
# Random effects model for multi-arm trials
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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

write.model(re_normal_gaus, "re-normal-gaus.txt")
MODELFILE.re <- c("re-normal-gaus.txt")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                          ---- Set up data for R2WinBUGs ----
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

nma_winbugs_datalist = function(data,treatments){

#Find column indexes
  t_loc = grep("t_|t.",colnames(data),fixed = F)
  y_loc = grep("y_|y.",colnames(data),fixed = F)
  se_loc = grep("se_|se.",colnames(data),fixed = F)
  na_loc = grep("na",colnames(data),fixed = F)
  v_loc = grep("V",colnames(data),fixed = F)
  
  
  t = as.matrix(pa_reac_wb[,t_loc])
  
  # number of treatments
  nt = length(treatments$description)
  
  y = as.matrix(cbind(rep(NA, length(t[,1])), data[,y_loc]))
  
  se = as.matrix(cbind(rep(NA, length(t[,1])), data[,se_loc]))
  
  na = as.vector(pa_reac_wb[,na_loc])
  
  V = as.vector(pa_reac_wb[,v_loc])
  
  ns2 = length(subset(na, na==2))
  ns3 = length(subset(na, na==3))
  ns4 = length(subset(na, na==4))
  
  data = list(nt=nt, ns2=ns2, ns3=ns3, ns4=ns4,t=t, y=y, se=se, na=na, V=V)
  
data
  
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                          ---- Run the analysis ----
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

nma_cont_fe = function(data,treatments){

  data = nma_winbugs_datalist(data,treatments)
  
  params.fe = c("meandif", 'SUCRA', 'best', 'totresdev', 'rk', 'dev', 'resdev', 'prob', "better")
  
  
  model = bugs(data, NULL, params.fe, model.file=MODELFILE.fe,
               n.chains = 3, n.iter = 60000, n.burnin = 40000, n.thin=1, 
               bugs.directory = "c:/Users/TheTimbot/Desktop/WinBUGS14", debug=F)

  comp = nma_cont_tables(model)
  
  sucra_rank = nma_cont_sucra_rank_prob(model)
#$$$$$$$$$$$$$$$$$$$$$$$$$
# All results
#$$$$$$$$$$$$$$$$$$$$$$$$$
  results = list(model,comp,rr,rankogram) 
  
  out = setNames(results, c("model", "Pairwise table","SUCRA and Ranking","Rankogram"))
  }


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                          ---- Create outputs ----
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

nma_cont_tables = function(model){
  
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  # Output all pairwise comparisons
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  

  
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

  comp
}
  

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



#                          ---- SUCRA, Rankings, and Probability Best Plots ----
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

nma_cont_sucra_rank_prob = function(model){
  
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


# Random-Effects: Probability summary
probs <- model$mean$prob; rownames(probs) <- treatments$description; colnames(probs) <- seq(1,nt)

dat <- melt(cbind(probs)); colnames(dat) <- c('Treatments', 'Rankings', 'Probability of Best Treatment')


rankogram = ggplot(dat,aes(x = Rankings, y = `Probability of Best Treatment`,fill = Treatments)) + 
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format())  +  ggtitle("Random Effects Rankogram") +scale_x_continuous(breaks=seq(1,nt))


  }

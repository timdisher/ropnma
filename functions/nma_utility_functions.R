#These functions are commonly used calculations for systematic reviews



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
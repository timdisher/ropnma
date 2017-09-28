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


#-----------------------------------------------------------------------------------------------------------------------------
# momlinc netgraph
#-----------------------------------------------------------------------------------------------------------------------------
# Takes a netmeta object and the int_char output from above to create a netgraph
# netmeta is a netmeta object e.g. pa_reac = netmeta(TE,seTE,treat1,treat2,studlab,data = some output from pairwise, sm = "MD")
# int_char is the table output from netmeta_xl_chars function
# outcome is a string that will be used for making names of files
# pointsize is the default pointsize for the graph
# location is a string specifying where results should be saved
#
# Example:
# pa_reac_contrast = pairwise(data = pa_reac,treat = trt.group, n= sample.size, mean = mean,sd = std.dev,studlab = study) 
# pa_reac_netconnect = netconnection(treat1,treat2,data = pa_reac_contrast)
# pa_reac_int = netmeta(TE,seTE,treat1,treat2,studlab,data = pa_reac_contrast, sm = "MD") ###required to drawn netgraph
# momlinc_netgraph(pa_reac_int,pa_reac_int_char,"pa_reac",2,location = "./figs/primary outcome/pain scales reactivity")


momlinc_netgraph = function(netmeta, int_char,pointsize){ ##outcome is a string
  
  
trt_details = tibble(Treatment = names(netmeta$TE.fixed[,1])) %>% left_join(int_char,by = "Treatment") %>% 
  mutate(w = `Number of Patients`/sum(`Number of Patients`)) ### normalizes weights based on sample size in treatment node


 netgraph(netmeta, col = rgb(116,104,170,maxColorValue = 225),
         points=TRUE, col.points = "aquamarine4", col.multiarm = "pink",cex.points=pointsize+trt_details$w*25,  ###weights point size by baseline + weight
         number.of.studies = FALSE, plastic = FALSE, thickness = "number.of.studies" ,
         multiarm = FALSE,
         lwd.max = 12,
         offset = 0.0575)


}

#-----------------------------------------------------------------------------------------------------------------------------
# Prepare data - may need to be tweaked for every data set until best approach finalized
#-----------------------------------------------------------------------------------------------------------------------------
#Example
#data = data_arm
#outcome = "PIPP
#timepoint.group = "reactivity"
#ref = "drops"

netmeta_prep = function(data,outcome,timepoint,ref){

data = data %>% filter(outcome == outcome, timepoint.group == timepoint) ##choose outcome

### converts data to correct format to allow for assessment of connectivity
pa_reac_contrast = pairwise(data = pa_reac,treat = treatment, n= sample.size, mean = mean,sd = std.dev,studlab = study) 

#- Assess whether network is connected -#
pa_reac_netconnect = netconnection(treat1,treat2,data = pa_reac_contrast)
pa_reac_int = netmeta(TE,seTE,treat1,treat2,studlab,data = pa_reac_contrast, sm = "MD")
}


#-----------------------------------------------------------------------------------------------------------------------------
# Placebo response graph- bar
#-----------------------------------------------------------------------------------------------------------------------------
#data = pa_reac
#ref = "drops"
#outcome = mean_imp


plac_resp_graph = function(data,x = "studlab",y = "mean",ref,outcome = "outcome",fill = "actual_timepoint",facet = FALSE,fv = "speculum",
                         ylab = "PIPP (mean or median)", type = "bar", size = 10, filter = "trt_group",
                         nrow = 2){

  if(type== "bar"){

  temp = data %>% filter(data[[filter]] == ref) %>% ggplot(aes_string(x = x, y = y,fill = fill)) + 
    geom_col() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(size = 13),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "PIPP Reactivity Drops Response Rate",
         x = "Study ID",
         y = ylab)+
    scale_fill_discrete(name = "Actual Timepoint")
  
  if(facet == FALSE){
   return(temp)}
  
  if(facet == TRUE){
    return(temp + facet_wrap(as.formula(paste("~",fv,sep="")), nrow = nrow) + labs(subtitle = paste("Faceted by",fv)))
    
  }
  }
  
  
  if(type == "scat"){

    temp = data %>% filter(data[[filter]] == ref) %>% ggplot(aes_string(x = x, y = y,colour = fill,size = size)) + 
      geom_point() +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            text = element_text(size = 13),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "PIPP Reactivity Drops Response Rate",
           x = "Study ID",
           y = ylab)+
      scale_fill_discrete(name = "Actual Timepoint") + geom_text(aes(label = studlab, size = 20),hjust =-0.15, vjust = 0)
    
    if(facet == FALSE){
return(temp)
}
    
    if(facet == TRUE){
return(temp + facet_wrap(as.formula(paste("~",fv,sep="")), nrow = nrow) + labs(subtitle = paste("Faceted by",fv)))
      

    }
  }
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

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# This function runs all pairwise meta analysis for an outcome
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# direct_comp_char = table output from net_metaxl tables
# data = contrast data output from the netmeta pairwise function
# list = a NULL list that will be populated with results
# outcome = string for titles and file names
# location = where it should be saved

# Example
# all_pairwise(pa_reac_direct_comp_chat,pa_reac_contrast,pa_reac_pairwise, outcome = "PIPP reactivity",location = getwd())

all_pairwise = function(direct_comp_char,data,location = getwd(),outcome,sm, cont = TRUE){

comps_data = direct_comp_char %>% filter(`# Studies`>1)
list = NULL
comp_matrix = str_split_fixed(comps_data$Treatment, " vs ",2)

if(cont == TRUE){
for(i in seq_along(comps_data$Treatment)){
  temp_data = data %>% filter(treat1 == comp_matrix[[i,2]] & treat2 == comp_matrix[i,1])
  list[[i]] = metacont(n2,mean2,sd2,n1,mean1,sd1, sm = sm, data = temp_data,studlab = studlab)
  names(list)[i] = comps_data$Treatment[i]
}
  } else{
  
  for(i in seq_along(comps_data$Treatment)){
    temp_data = data %>% filter(treat1 == comp_matrix[[i,2]] & treat2 == comp_matrix[i,1])
    list[[i]] = metabin(n2,event2,sd2,n1,event1,sd1, sm = sm, data = temp_data,studlab = studlab)
    names(list)[i] = comps_data$Treatment[i] 
  
  
}
}


for(i in seq_along(list)){
  pdf(paste(location,"/",outcome,"_",names(list[i]),"_pw_ma.pdf",sep=""),width = 12, height = 6)
  forest(list[[i]], lab.e = comp_matrix[i,1], lab.c = comp_matrix[i,2])
  grid.text(paste(outcome,names(list[i]),sep = " "), .5, .8, gp=gpar(cex=2))
  dev.off()
}

list
}


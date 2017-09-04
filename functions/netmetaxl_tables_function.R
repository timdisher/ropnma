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

netmeta_xl_chars = function(data,outcome,ref,treat = "treatment",location = getwd()){

dir.create(location,recursive = TRUE)
  
contrast = pairwise(data = data,treat = data[[treat]], n= sample.size, mean = mean,sd = std.dev,studlab = study)

studies = (data %>% select(study) %>% distinct() %>% count())[[1]] ## count the number of studies
trts = (data %>% select(treat) %>% distinct() %>% count())[[1]] ## count the number of treatments
totn = data %>% select(sample.size) %>% sum() # count the number of patients
poss_pw = choose(trts,2) # total possible pairwise comparisons


#-- Create a numbered list of treatments with ref as #1--#
data[[treat]] = fct_relevel(data[[treat]], ref)
names= as.character((data %>% select(treat) %>% distinct() %>% arrange_("trt.group"))[[treat]])

names = tibble("Treatment Number" = seq(1:trts),
               "Treatment Description" = names)

               
#- Assign treatment numbers based on tibble above -#
contrast = names %>% rename(treat2_asnum = `Treatment Number`,
                            treat2 = `Treatment Description`) %>% right_join(contrast, by = "treat2")
contrast = names %>% rename(treat1_asnum = `Treatment Number`,
                            treat1 = `Treatment Description`) %>% right_join(contrast, by = "treat1")


assign(paste(outcome,"_data_nma",sep=""), contrast,envir=globalenv())
               

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

assign(paste(outcome,"_net_char",sep=""), net_char,envir=globalenv())
stargazer(net_char, summary = FALSE,title = "Network Characteristics",
          out = paste(location,"/",outcome,"_net_char.html",sep=""),
          rownames = FALSE)
#-----------------------------------
# Intervention characteristics table
#-----------------------------------,
style = "asq"
int_nums = data %>% group_by_(treat) %>% summarise("Number of Comparisons" = n(),
                                                   "Number of Patients" = sum(sample.size)) %>% rename_(Treatment = treat)


int_char <- tibble("Treatment" = names[[2]]) %>% left_join(int_nums)
  
assign(paste(outcome,"_int_char",sep=""), int_char,envir=globalenv())
stargazer(int_char, summary = FALSE,title = "Intervention Characteristics",
          out = paste(location,"/",outcome,"_int_char.html",sep=""),
          rownames = FALSE)

#-----------------------------------
# Direct comparison characteristics table
#-----------------------------------
direct_comp_char <- direct %>% filter(nstud > 0) %>% left_join(names, by = c("trt1" = "Treatment Number")) %>%
left_join(names, by = c("trt2" = "Treatment Number")) %>% rename(treatment_1 = `Treatment Description.x`,
                                                                 treatment_2 = `Treatment Description.y`) %>%
  mutate(Treatment = paste(treatment_1,"vs",treatment_2)) %>%
  select(Treatment,nstud,ntot) %>% rename("# Studies" = nstud,
                                          "# Patients" = ntot)

assign(paste(outcome,"_direct_comp_char",sep=""), direct_comp_char,envir=globalenv())
stargazer(direct_comp_char, summary = FALSE,title = "Direct Comparison Characteristics",
          out = paste(location,"/",outcome,"_direct_comp_char.html",sep=""),
          rownames = FALSE)


print(net_char)
print(int_char)
print(direct_comp_char)
}


#-----------------------------------------------------------------------------------------------------------------------------
# momlinc netgraph
#-----------------------------------------------------------------------------------------------------------------------------
# Takes a netmeta object and the int_char output from above to create a netgraph
#
#
momlinc_netgraph = function(netmeta, int_char,outcome,pointsize,location = getwd()){ ##outcome is a string
  
  dir.create(location,recursive = TRUE)
  
  
trt_details = tibble(Treatment = names(netmeta$TE.fixed[,1])) %>% left_join(int_char,by = "Treatment") %>% 
  mutate(w = `Number of Patients`/sum(`Number of Patients`)) ### normalizes weights based on sample size in treatment node

temp_graph = netgraph(netmeta, col = rgb(116,104,170,maxColorValue = 225),
         points=TRUE, col.points = "aquamarine4", col.multiarm = "pink",cex.points=pointsize+trt_details$w*25,  ###weights point size by baseline + weight
         number.of.studies = FALSE, plastic = FALSE, thickness = "number.of.studies" ,
         multiarm = FALSE,
         lwd.max = 12,
         offset = 0.0575)

pdf(paste(location,"/",outcome,"_netgraph.pdf",sep=""),width = 10)
netgraph(netmeta, col = rgb(116,104,170,maxColorValue = 225),
         points=TRUE, col.points = "aquamarine4", col.multiarm = "pink",cex.points=pointsize+trt_details$w*25,  ###weights point size by baseline + weight
         number.of.studies = FALSE, plastic = FALSE, thickness = "number.of.studies" ,
         multiarm = FALSE,
         lwd.max = 12,
         offset = 0.0575)
dev.off()
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


plac_resp_graph = function(data,x = "study",y = "mean_imp",ref,outcome = "outcome",fill = "actual.timepoint",facet = FALSE,fv = "speculum",
                         ylab = "PIPP (mean or median)", type = "bar", size = 10, filter = "trt.group",
                         location = getwd(),
                         width = 7,
                         height  = 7,
                         nrow = 2){

  if(type== "bar"){
  dir.create(location,recursive = TRUE)
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
    print(temp)
    ggsave(paste(location,"/",outcome,"_plac_resp_bar.pdf",sep=""),device = "pdf", height = height, width = width)}
  
  if(facet == TRUE){
    print(temp + facet_wrap(as.formula(paste("~",fv,sep="")), nrow = nrow) + labs(subtitle = paste("Faceted by",fv)))
    
    ggsave(paste(location,"/",outcome,"_plac_resp_scat_facet_",fv,".pdf",sep=""),device = "pdf",height = height, width = width)
  }
  }
  
  
  if(type == "scat"){
    dir.create(location,recursive = TRUE)
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
      print(temp)
      ggsave(paste(location,"/",outcome,"_plac_resp_scat.pdf",sep=""),device = "pdf", height = height, width = width)}
    
    if(facet == TRUE){
      print(temp + facet_wrap(as.formula(paste("~",fv,sep="")), nrow = nrow) + labs(subtitle = paste("Faceted by",fv)))
      
      ggsave(paste(location,"/",outcome,"_plac_resp_scat_facet_",fv,".pdf",sep=""),device = "pdf",height = height, width = width)
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




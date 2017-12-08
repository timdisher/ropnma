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
##requires tidyverse
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



#============================================================================== =
# Pairwise function from the netmeta package (copied here to avoid loading ====== 
# package for on function/protect against changes in future versions) ========= =
# ============================================================================= =

pairwise = function (treat, event, n, mean, sd, TE, seTE, time, data = NULL, 
                     studlab, incr = 0.5, allincr = FALSE, addincr = FALSE, allstudies = FALSE, 
                     ...) 
{
  if (is.null(data)) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  studlab <- eval(mf[[match("studlab", names(mf))]], data, 
                  enclos = sys.frame(sys.parent()))
  treat <- eval(mf[[match("treat", names(mf))]], data, enclos = sys.frame(sys.parent()))
  event <- eval(mf[[match("event", names(mf))]], data, enclos = sys.frame(sys.parent()))
  n <- eval(mf[[match("n", names(mf))]], data, enclos = sys.frame(sys.parent()))
  mean <- eval(mf[[match("mean", names(mf))]], data, enclos = sys.frame(sys.parent()))
  sd <- eval(mf[[match("sd", names(mf))]], data, enclos = sys.frame(sys.parent()))
  TE <- eval(mf[[match("TE", names(mf))]], data, enclos = sys.frame(sys.parent()))
  seTE <- eval(mf[[match("seTE", names(mf))]], data, enclos = sys.frame(sys.parent()))
  time <- eval(mf[[match("time", names(mf))]], data, enclos = sys.frame(sys.parent()))
  args <- list(...)
  nam.args <- names(args)
  if (is.null(treat)) 
    stop("Argument 'treat' mandatory.")
  if (is.list(treat)) 
    chklist(treat)
  if (!is.null(event)) 
    if (is.list(event)) 
      chklist(event)
  else meta:::chknumeric(event)
  if (!is.null(n)) 
    if (is.list(n)) 
      chklist(n)
  else meta:::chknumeric(n)
  if (!is.null(mean)) 
    if (is.list(mean)) 
      chklist(mean)
  else meta:::chknumeric(mean)
  if (!is.null(sd)) 
    if (is.list(sd)) 
      chklist(sd)
  else meta:::chknumeric(sd)
  if (!is.null(TE)) 
    if (is.list(TE)) 
      chklist(TE)
  else meta:::chknumeric(TE)
  if (!is.null(seTE)) 
    if (is.list(seTE)) 
      chklist(seTE)
  else meta:::chknumeric(seTE)
  if (!is.null(time)) 
    if (is.list(time)) 
      chklist(time)
  else meta:::chknumeric(time)
  meta:::chknumeric(incr, min = 0, single = TRUE)
  meta:::chklogical(allincr)
  meta:::chklogical(addincr)
  meta:::chklogical(allstudies)
  if (!is.null(event) & !is.null(n) & is.null(mean) & is.null(sd) & 
      is.null(TE) & is.null(seTE) & is.null(time)) 
    type <- "binary"
  else if (is.null(event) & !is.null(n) & !is.null(mean) & 
           !is.null(sd) & is.null(TE) & is.null(seTE) & is.null(time)) 
    type <- "continuous"
  else if (!is.null(event) & is.null(n) & is.null(mean) & 
           is.null(sd) & is.null(TE) & is.null(seTE) & !is.null(time)) 
    type <- "count"
  else if (is.null(event) & is.null(n) & is.null(mean) & is.null(sd) & 
           !is.null(TE) & !is.null(seTE) & is.null(time)) 
    type <- "generic"
  else stop("Type of outcome unclear. Please provide the necessary information:\n  - event, n (binary outcome)\n  - n, mean, sd (continuous outcome)\n  - TE, seTE (generic outcome)\n  - event, time (incidence rates).")
  treat.list <- list()
  event.list <- list()
  n.list <- list()
  mean.list <- list()
  sd.list <- list()
  TE.list <- list()
  seTE.list <- list()
  time.list <- list()
  if (type == "binary") {
    listformat <- is.list(event) & is.list(n)
    if (!listformat) {
      if (is.null(studlab)) 
        stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      tdat <- data.frame(studlab, treat, event, n, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", 
                        all.x = TRUE)
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        n.list[[i]] <- tres.i$n
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      treat <- treat.list
      event <- event.list
      n <- n.list
    }
  }
  else if (type == "continuous") {
    listformat <- is.list(n) & is.list(mean) & is.list(sd)
    if (!listformat) {
      if (is.null(studlab)) 
        stop("Argument 'studlab' mandatory if argument 'mean' is a vector.")
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      tdat <- data.frame(studlab, treat, n, mean, sd, 
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", 
                        all.x = TRUE)
        treat.list[[i]] <- tres.i$treat
        n.list[[i]] <- tres.i$n
        mean.list[[i]] <- tres.i$mean
        sd.list[[i]] <- tres.i$sd
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      treat <- treat.list
      n <- n.list
      mean <- mean.list
      sd <- sd.list
    }
  }
  else if (type == "count") {
    listformat <- is.list(event) & is.list(time)
    if (!listformat) {
      if (is.null(studlab)) 
        stop("Argument 'studlab' mandatory if argument 'event' is a vector.")
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      tdat <- data.frame(studlab, treat, event, time, 
                         stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", 
                        all.x = TRUE)
        treat.list[[i]] <- tres.i$treat
        event.list[[i]] <- tres.i$event
        time.list[[i]] <- tres.i$time
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      treat <- treat.list
      event <- event.list
      time <- time.list
    }
  }
  else if (type == "generic") {
    listformat <- is.list(TE) & is.list(seTE)
    if (!listformat) {
      if (is.null(studlab)) 
        stop("Argument 'studlab' mandatory if argument 'TE' is a vector.")
      ttab <- table(as.character(studlab), as.character(treat))
      n.arms <- apply(ttab, 1, sum)
      tdat <- data.frame(studlab, treat, TE, seTE, stringsAsFactors = FALSE)
      tdat <- tdat[order(tdat$studlab, tdat$treat), ]
      studlab <- names(n.arms)
      tres <- data.frame(studlab = studlab, stringsAsFactors = FALSE)
      for (i in 1:max(n.arms)) {
        tdat.i <- tdat[!duplicated(tdat$studlab), ]
        tres.i <- merge(tres, tdat.i, by = "studlab", 
                        all.x = TRUE)
        treat.list[[i]] <- tres.i$treat
        TE.list[[i]] <- tres.i$TE
        seTE.list[[i]] <- tres.i$seTE
        tdat <- tdat[duplicated(tdat$studlab), ]
      }
      treat <- treat.list
      TE <- TE.list
      seTE <- seTE.list
    }
  }
  if (is.null(studlab)) 
    studlab <- seq(along = treat[[1]])
  if (length(studlab) != length(unique(studlab))) 
    stop("Study labels must all be distinct.")
  levs <- studlab
  narms <- length(treat)
  nstud <- length(studlab)
  sumzero <- function(x) sum(x[!is.na(x)] == 0)
  if (type == "binary") {
    if (length(event) != narms) 
      stop("Different length of lists 'treat' and 'event'.")
    if (length(n) != narms) 
      stop("Different length of lists 'treat' and 'n'.")
    n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 
                     1, sumzero)
    n.all <- apply(matrix(unlist(n), ncol = length(event)) - 
                     matrix(unlist(event), ncol = length(event)), 1, 
                   sumzero)
    incr.study <- rep(0, length(n.zeros))
    if ("sm" %in% nam.args) 
      sm <- args$sm
    else sm <- gs("smbin")
    sm <- meta:::setchar(sm, c("OR", "RD", "RR", "ASD"))
    sparse <- switch(sm, OR = (n.zeros > 0) | (n.all > 0), 
                     RD = (n.zeros > 0) | (n.all > 0), RR = (n.zeros > 
                                                               0) | (n.all > 0), ASD = rep(FALSE, length(n.zeros)))
    if (!allincr & !addincr) 
      incr.study[sparse] <- incr
    else if (addincr) 
      incr.study[] <- incr
    else {
      if (any(n.zeros > 0)) 
        incr.study[] <- incr
      else incr.study[] <- 0
    }
    for (i in 1:(narms - 1)) {
      if (i == 1 & (length(treat[[i]]) != length(event[[i]]))) 
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i == 1 & (length(event[[i]]) != length(n[[i]]))) 
        stop("Different length of element ", i, " of lists 'event' and 'n'.")
      for (j in (i + 1):narms) {
        if (length(treat[[j]]) != length(event[[j]])) 
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(event[[j]]) != length(n[[j]])) 
          stop("Different length of element ", j, " of lists 'event' and 'n'.")
        dat <- data.frame(TE = NA, seTE = NA, studlab = studlab, 
                          treat1 = treat[[i]], treat2 = treat[[j]], 
                          event1 = event[[i]], n1 = n[[i]], event2 = event[[j]], 
                          n2 = n[[j]], incr = incr.study, allstudies = allstudies, 
                          stringsAsFactors = FALSE)
        dat <- dat[!(is.na(dat$event1) & is.na(dat$n1)), 
                   ]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$n2)), 
                   ]
        if (nrow(dat) > 0) {
          m1 <- metabin(dat$event1, dat$n1, dat$event2, 
                        dat$n2, incr = dat$incr, addincr = TRUE, 
                        allstudies = allstudies, ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | 
                           dat$seTE <= 0, ]
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else if (i == 1 & j == 2) 
          stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  if (type == "continuous") {
    if (length(n) != narms) 
      stop("Different length of lists 'treat' and 'n'.")
    if (length(mean) != narms) 
      stop("Different length of lists 'treat' and 'mean'.")
    if (length(sd) != narms) 
      stop("Different length of lists 'treat' and 'sd'.")
    for (i in seq_len(narms)) {
      if (length(treat[[i]]) != length(n[[i]])) 
        stop("Different length of element ", i, " of lists 'treat' and 'n'.", 
             call. = FALSE)
      if (length(treat[[i]]) != length(mean[[i]])) 
        stop("Different length of element ", i, " of lists 'treat' and 'mean'.", 
             call. = FALSE)
      if (length(treat[[i]]) != length(sd[[i]])) 
        stop("Different length of element ", i, " of lists 'treat' and 'sd'.", 
             call. = FALSE)
      if (length(treat[[i]]) != nstud) 
        stop("Different length of study labels and element ", 
             i, " of list 'treat'.", call. = FALSE)
    }
    if ("sm" %in% nam.args && (tolower(args$sm) == "smd" & 
                               narms > 2)) {
      pooled.sd <- function(sd, n) {
        sel <- !is.na(sd) & !is.na(n)
        if (any(sel)) 
          res <- sqrt(sum((n[sel] - 1) * sd[sel]^2)/sum(n[sel] - 
                                                          1))
        else res <- NA
        res
      }
      N <- matrix(unlist(n), ncol = narms, nrow = nstud, 
                  byrow = FALSE)
      M <- matrix(unlist(mean), ncol = narms, nrow = nstud, 
                  byrow = FALSE)
      S <- matrix(unlist(sd), ncol = narms, nrow = nstud, 
                  byrow = FALSE)
      sel.n <- apply(!is.na(N) & N > 0, 1, sum) > 2
      sel.mean <- apply(!is.na(M), 1, sum) > 2
      sel.sd <- apply(!is.na(S) & S > 0, 1, sum) > 2
      sel <- sel.n & sel.mean & sel.sd
      if (any(sel)) {
        N <- N[sel, , drop = FALSE]
        S <- S[sel, , drop = FALSE]
        sd.p <- rep_len(NA, nrow(N))
        for (i in seq_len(nrow(N))) sd.p[i] <- pooled.sd(S[i, 
                                                           ], N[i, ])
      }
      for (i in seq_len(narms)) sd[[i]][sel] <- ifelse(is.na(sd[[i]][sel]), 
                                                       NA, sd.p)
    }
    for (i in 1:(narms - 1)) {
      for (j in (i + 1):narms) {
        dat <- data.frame(TE = NA, seTE = NA, studlab = studlab, 
                          treat1 = treat[[i]], treat2 = treat[[j]], 
                          n1 = n[[i]], mean1 = mean[[i]], sd1 = sd[[i]], 
                          n2 = n[[j]], mean2 = mean[[j]], sd2 = sd[[j]], 
                          stringsAsFactors = FALSE)
        dat <- dat[!(is.na(dat$n1) & is.na(dat$mean1) & 
                       is.na(dat$sd1)), ]
        dat <- dat[!(is.na(dat$n2) & is.na(dat$mean2) & 
                       is.na(dat$sd2)), ]
        if (nrow(dat) > 0) {
          m1 <- metacont(dat$n1, dat$mean1, dat$sd1, 
                         dat$n2, dat$mean2, dat$sd2, ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | 
                           dat$seTE <= 0, ]
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else if (i == 1 & j == 2) 
          stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  if (type == "generic") {
    if (length(TE) != narms) 
      stop("Different length of lists 'treat' and 'TE'.")
    if (length(seTE) != narms) 
      stop("Different length of lists 'treat' and 'seTE'.")
    for (i in 1:(narms - 1)) {
      if (i == 1 & (length(treat[[i]]) != length(TE[[i]]))) 
        stop("Different length of element ", i, " of lists 'treat' and 'TE'.")
      if (i == 1 & (length(treat[[i]]) != length(seTE[[i]]))) 
        stop("Different length of element ", i, " of lists 'treat' and 'seTE'.")
      for (j in (i + 1):narms) {
        if (length(treat[[j]]) != length(TE[[j]])) 
          stop("Different length of element ", j, " of lists 'treat' and 'TE'.")
        if (length(treat[[j]]) != length(seTE[[j]])) 
          stop("Different length of element ", j, " of lists 'treat' and 'seTE'.")
        dat <- data.frame(TE = NA, seTE = NA, studlab = studlab, 
                          treat1 = treat[[i]], treat2 = treat[[j]], 
                          TE1 = TE[[i]], seTE1 = seTE[[i]], TE2 = TE[[j]], 
                          seTE2 = seTE[[j]], stringsAsFactors = FALSE)
        dat <- dat[!(is.na(dat$TE1) & is.na(dat$seTE1)), 
                   ]
        dat <- dat[!(is.na(dat$TE2) & is.na(dat$seTE2)), 
                   ]
        if (nrow(dat) > 0) {
          m1 <- metagen(dat$TE1 - dat$TE2, sqrt(dat$seTE1^2 + 
                                                  dat$seTE2^2), ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | 
                           dat$seTE <= 0, ]
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else if (i == 1 & j == 2) 
          stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  if (type == "count") {
    if (length(event) != narms) 
      stop("Different length of lists 'treat' and 'event'.")
    if (length(time) != narms) 
      stop("Different length of lists 'treat' and 'time'.")
    n.zeros <- apply(matrix(unlist(event), ncol = length(event)), 
                     1, sumzero)
    incr.study <- rep(0, length(n.zeros))
    sparse <- n.zeros > 0
    if (!allincr & !addincr) 
      incr.study[sparse] <- incr
    else if (addincr) 
      incr.study[] <- incr
    else {
      if (any(n.zeros > 0)) 
        incr.study[] <- incr
      else incr.study[] <- 0
    }
    for (i in 1:(narms - 1)) {
      if (i == 1 & (length(treat[[i]]) != length(event[[i]]))) 
        stop("Different length of element ", i, " of lists 'treat' and 'event'.")
      if (i == 1 & (length(treat[[i]]) != length(time[[i]]))) 
        stop("Different length of element ", i, " of lists 'treat' and 'time'.")
      for (j in (i + 1):narms) {
        if (length(treat[[j]]) != length(event[[j]])) 
          stop("Different length of element ", j, " of lists 'treat' and 'event'.")
        if (length(treat[[j]]) != length(time[[j]])) 
          stop("Different length of element ", j, " of lists 'treat' and 'time'.")
        dat <- data.frame(TE = NA, seTE = NA, studlab = studlab, 
                          treat1 = treat[[i]], treat2 = treat[[j]], 
                          event1 = event[[i]], time1 = time[[i]], event2 = event[[j]], 
                          time2 = time[[j]], incr = incr.study, stringsAsFactors = FALSE)
        dat <- dat[!(is.na(dat$event1) & is.na(dat$time1)), 
                   ]
        dat <- dat[!(is.na(dat$event2) & is.na(dat$time2)), 
                   ]
        if (nrow(dat) > 0) {
          m1 <- metainc(dat$event1, dat$time1, dat$event2, 
                        dat$time2, incr = dat$incr, addincr = TRUE, 
                        allstudies = allstudies, ...)
          dat$TE <- m1$TE
          dat$seTE <- m1$seTE
          dat.NAs <- dat[is.na(dat$TE) | is.na(dat$seTE) | 
                           dat$seTE <= 0, ]
          if (i == 1 & j == 2) {
            res <- dat
            res.NAs <- dat.NAs
          }
          else {
            res <- rbind(res, dat)
            res.NAs <- rbind(res.NAs, dat.NAs)
          }
        }
        else if (i == 1 & j == 2) 
          stop("No studies available for comparison of first and second treatment.")
      }
    }
  }
  sel.treat <- as.character(res$treat1) == as.character(res$treat2)
  if (any(sel.treat)) {
    stop(paste("Identical treatments for the following studies:\n  ", 
               paste(paste("'", studlab[sel.treat], "'", sep = ""), 
                     collapse = " - "), sep = ""))
  }
  sel.study <- !(studlab %in% unique(as.character(res$studlab)))
  if (any(sel.study)) 
    warning(paste("The following studies are not considered in the analysis\n  ", 
                  "(due to single study arm or missing values):\n  ", 
                  paste(paste("'", studlab[sel.study], "'", sep = ""), 
                        collapse = " - "), sep = ""))
  if (nrow(res.NAs) > 0) {
    warning("Comparison", if (nrow(res.NAs) > 1) 
      "s", " with missing TE / seTE or zero seTE", " will not be considered in network meta-analysis.", 
      call. = FALSE)
    cat(paste("Comparison", if (nrow(res.NAs) > 1) 
      "s", " will not be considered in network meta-analysis:\n", 
      sep = ""))
    prmatrix(res.NAs, quote = FALSE, right = TRUE, na.print = "NA", 
             rowlab = rep("", nrow(res.NAs)))
  }
  attr(res, "sm") <- m1$sm
  attr(res, "method") <- m1$method
  attr(res, "version") <- packageDescription("netmeta")$Version
  res <- res[order(factor(res$studlab, levels = levs), res$treat1, 
                   res$treat2), ]
  rownames(res) <- 1:nrow(res)
  class(res) <- c(class(res), "pairwise")
  res
}
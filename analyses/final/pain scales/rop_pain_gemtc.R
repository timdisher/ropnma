#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Network Meta-Analysis in R - Adapted from script provided by Cornerstone Research Group

# Project: Pain from Retinopathy of Prematurity Eye Exams
# Outcome: Pain reactivity

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(forestplot)
library(personograph)

source("./analyses/final/pain scales/rop_explore_pipp.R")
source("./functions/gemtc test/nma_cont_gemtc.R")
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Load data in WInBugs Format
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# #========================================================================================
# 
# 
# Outcome: Pain Reactivity-----
# 
# #========================================================================================

#=============
# Functions for this project only
#==============

#Create and run gemtc analyses====

data = pa_recov_data$pa$gemtc$data
set_net = function(data = pa_reac_data$pa$gemtc$data, model = "consistency",regressor = NULL, type = "random"){
  network = mtc.network(data.re = data[,1:4],
                        studies = data[,c(1,5:8)])
  
  if(model == "consistency"){
  model = mtc.model(network, type = "consistency",
                    linearModel = type,likelihood = "normal",
                    link = "identity")
  } else{
    model = mtc.model(network, type = "regression",
                      linearModel = type,likelihood = "normal",
                      regressor = regressor,
                      link = "identity")
  }
  results = mtc.run(model)
  
  sucleague = suc_leag_rop(results)
  
  results = list(network = network,results = results, suc = sucleague$suc,league = sucleague$league)
  
  
}

# Get sucra and league tables=====

suc_leag_rop = function(data = pa_reac_data$pa$results){
  
  suc = as.data.frame(sucra(data, direction = -1)) %>% rownames_to_column("treat") %>% rename(sucra = `sucra(data, direction = -1)`) %>% 
    mutate(treat = paste("d.drops.",treat,sep = ""))%>% arrange(-sucra)
  
  basicp = relative.effect(data,t1 = c("drops"),preserve.extra = FALSE)
  basicp = as.data.frame(as.matrix(basicp$samples)) %>% mutate(d.drops.drops = 0)
  
  league = league(results = basicp, order = suc)
  
  list(sucra = suc, league = league)
}


# Get quick mcmc diagnostics=======
gemtc_diag = function(results = pa_reac_data$pa$results){
  windows(record = TRUE)
  plot(results)
  gelman.plot(results)
  gelman.diag(results)
}


# Conduct power analysis based on posterior SD and mid

power_table = function(data = pa_reac_data$sa6$results, direct_comps = chars$direct_zeros, mid = 2){
  power = direct_comps %>% rename(ctrl = `Treatment Description.x`,
                                  trt = `Treatment Description.y`) %>%  unite(comp,trt,ctrl,sep = " vs ") %>% select(comp, ntot,nstud) %>% mutate(num = seq(1,length(comp),1))
  
  retrodesign <- function(A, s, alpha=.05, df=Inf, n.sims=10000) { 
    z <- qt(1-alpha/2, df)
    p.hi <- 1 - pt(z-A/s, df)
    p.lo <- pt(-z-A/s, df)  
    power <- p.hi + p.lo
    typeS <- p.lo/power
    estimate <- A + s*rt(n.sims,df)
    significant <- abs(estimate) > s*z
    exaggeration <- mean(abs(estimate)[significant])/A
    return(list(power=power,
                typeS=typeS, exaggeration=exaggeration))
  }
  
  outcome_names = power[grep("vs drops$",power$comp),c(1)] %>% separate(comp,c("t2","t1"), sep = " vs ") %>% select(t2)
  
  pa_reac_pt_res = data.frame(outcome = outcome_names$t2,
                              median = 0,
                              low = 0,
                              high = 0,
                              sd = 0)
  for(i in seq_along(outcome_names$t2)){
    temp = summary(relative.effect(data, t1 = "drops",t2 = outcome_names$t2[i])) 
    
    pa_reac_pt_res$median[i] = temp$summaries$quantiles[1,3]
    pa_reac_pt_res$low[i] = temp$summaries$quantiles[1,1]
    pa_reac_pt_res$high[i] = temp$summaries$quantiles[1,5]
    pa_reac_pt_res$sd[i] = temp$summaries$statistics[1,2]
    
    
  }
  
  
  power_table = pa_reac_pt_res %>% mutate(gelman_power = 0)
  
  for(i in seq_along(power_table$outcome)){
    power_table$gelman_power[i] = retrodesign(mid,power_table$sd[i])[["power"]]
    
  }
  power_table
}

# Create forest plots including power======
momlinc_fp_p = function(data = pa_reac_power, names = comp_names, width = 7, height = 5){
  
  pa_reac_forest_data = data %>% mutate(outcome = names) %>% arrange(median) 
  
  pa_reac_plot_data = pa_reac_forest_data %>% select(median, low, high) %>% rename(mean = median,
                                                                                   lower = low,
                                                                                   upper = high)
  pa_reac_table_data = pa_reac_forest_data %>%  mutate(mean_cri = paste(round(median,2)," (",round(low,2)," to ",round(high,2),")",sep=""),
                                                       gelman_power = round(pa_reac_forest_data$gelman_power,2)) %>% select(-median,-low,-high,-sd)
  
  
  windows(width = width, height = height)
  forestplot(rbind(c("Comparison","Power","Mean Difference (95% CrI)"),pa_reac_table_data),
             c(NA,as.numeric(as.character(pa_reac_plot_data$mean))),
             c(NA,as.numeric(as.character(pa_reac_plot_data$lower))),
             c(NA,as.numeric(as.character(pa_reac_plot_data$upper))),
             graph.pos = 3, graphwidth = unit(50,'mm'),
             is.summary = c(TRUE,rep(FALSE,12)),
             hrzl_lines = gpar(col="#444444"),
             align = c("l",rep("c",5),"l"),
             boxsize = 0.4,
             colgap = unit(3,"mm"),
             lineheight = unit(7,"mm"),
             txt_gp = fpTxtGp(label = gpar(fontsize = 11, family = "calibri")),
             col = fpColors(box = "mediumpurple", line = "midnightblue"),
             xlab = "Compared to drops alone",
             title = "PIPP Reactivity (Intervention vs Anesthetic eye drops alone)"
  )
  
}

# Create heatplots for ROP outcomes====

rop_heatplots = function(data = pa_reac_data,analyses_index = analyses,pub_names = names, treats = treat_names){
  
  suc_list = NULL
  
  for(i in 1:length(analyses_index)){
    suc_list[[i]] = pa_reac_data[[analyses_index[[i]]]]$mod$suc
  }
  
  sa_table = suc_list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="treat"), .)
  
  colnames(sa_table) = c("Treatment",names)
  
  sa_table["Treatment"] = treat_names
  
  
  plot = heatplot(sa_table)
  
  list(data = suc_list, plot = plot)
  
}

calc_n = function(data = pa_reac_data$pa$gemtc$data,input = pa_reac_data$pa$gemtc$input){
  
inc = data$study %>% unique()

n = input %>% filter(studlab %in% inc) %>% group_by(studlab) %>% summarize(n = sum(n))
studies = length(n$n)
n = sum(n$n)
list(studies = studies, n = n)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Prep Data
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_reac_data = NULL

pa_reac_data$pa$gemtc = prep_gem(pa_reac)

pa_reac_data$pa$gemtc$data = pa_reac_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","avg_pma","pub_type","oa_rob_sub")],by = c("study" = "studlab"))%>% 
  left_join(pa_reac[c("studlab","imputed_mean")] %>% distinct(), by = c("study" = "studlab")) %>%
  mutate(pub_type = ifelse(pub_type == "journal",0,1),
         oa_rob_sub = ifelse(oa_rob_sub == "low",0,1),
         imputed_mean = ifelse(imputed_mean == "no",0,1)) %>% replace_na(list(imputed_mean = 0))


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




pa_reac_data$pa$mod = set_net()

gemtc_diag(pa_reac_data$pa$mod$results)

summary(pa_reac_data$pa$mod$results)


#Analysis of heterogeneity and rough plot of effects versus reference 
# Nodeplit (primary)
pa_reac_nodesplit = mtc.nodesplit(pa_reac_data$pa$mod$network)
as.data.frame(summary(pa_reac_nodesplit))

plot(summary(pa_reac_nodesplit))
#ANOHE - Fits a UME, USE, and Cons model... good visualization
pa_reac_data$pa$anohe = mtc.anohe(pa_reac_data$pa$mod$network)

pdf("./figs/pain scales reactivity/pa_reac_anohe.pdf", height = 11, width = 8.5)
plot(summary(pa_reac_data$pa$anohe))
dev.off()


#Write league table to file
write.csv(pa_reac_data$pa$mod$league,"./tables/pa_reac_league.csv")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- Drop Boyle
# Rationale: Boyle is the source of large amounts of heterogeneity and the trial is at high risk of bias + has extremely small sample per arm'
# Cannot be assessed through meta-regression since it is the only study in the inconsistent comparison. This dataset used for all following analyses.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#Drop Boyle
pa_reac_data$pa$gemtc$data_nb =  pa_reac_data$pa$gemtc$data %>% filter(study != "Boyle 2006")

pa_reac_data$sa1$mod =  set_net(pa_reac_data$pa$gemtc$data_nb)

gemtc_diag(pa_reac_data$sa1$mod$results)

summary(pa_reac_data$sa1$mod$results)

#Nodesplit
pa_reac_data$sa1$nodesplit = mtc.nodesplit(pa_reac_data$sa1$mod$network)

summary(pa_reac_data$sa1$nodesplit)

#ANOHE
pa_reac_data$sa1$anohe = mtc.anohe(pa_reac_data$sa1$mod$network)
plot(summary(pa_reac_data$sa1$anohe))

summary(pa_reac_data$sa1$mod$results)

gemtc::forest(relative.effect.table(pa_reac_data$sa1$mod$results,"drops"))


# Add Boyle to list of excluded studies
boyle_n = pa_reac_data$pa$gemtc$input %>% filter(studlab == "Boyle 2006") %>%
  summarise(sample = sum(n),
            treat = paste(treatment,collapse=' vs ')) %>%
  mutate(reason = "High risk of bias, source of significant inconsistency")

pa_reac_excluded = rbind(pa_reac_excluded,boyle_n)
  
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2
# ----- Meta-regression on pma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sa2_regressor = list(coefficient = "shared",
                     variable = "avg_pma",
                     control = "drops")

pa_reac_data$sa2$mod = set_net(data = pa_reac_data$pa$gemtc$data_nb %>% drop_na(avg_pma) %>% droplevels("treatment"),
        model = "regression",regressor = sa2_regressor)

gemtc_diag(pa_reac_data$sa2$mod$results)

summary(pa_reac_data$sa2$mod$results)

forest(relative.effect.table(pa_reac_data$sa2$mod$results),"drops")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3
# ----- MR Risk of bias
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sa3_regressor = list(coefficient = "shared",
                     variable = "oa_rob_sub",
                     control = "drops")

pa_reac_data$sa3$mod = set_net(data = pa_reac_data$pa$gemtc$data_nb,
                               model = "regression",regressor = sa3_regressor)

gemtc_diag(pa_reac_data$sa3$mod$results)

summary(pa_reac_data$sa3$mod$results)

forest(relative.effect.table(pa_reac_data$sa3$mod$results),"drops")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 4
# ----- Meta-regression on pub type, not enough data, no convergence, just removing instead
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_reac_data$sa4$mod = set_net(pa_reac_data$pa$gemtc$data_nb %>% filter(pub_type == 0))

calc_n(data = pa_reac_data$pa$gemtc$data_nb %>% filter(pub_type == 0),input = pa_reac_data$pa$gemtc$input)

gemtc_diag(pa_reac_data$sa4$mod$results)


summary(pa_reac_data$sa4$mod$results)

forest(relative.effect.table(pa_reac_data$sa4$mod$results),"drops")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 5
# ----- Remove imputed mean
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_reac_data$sa5$mod = set_net(pa_reac_data$pa$gemtc$data_nb %>% filter(imputed_mean == 0))

calc_n(data = pa_reac_data$pa$gemtc$data_nb %>% filter(imputed_mean == 0),input = pa_reac_data$pa$gemtc$input)


gemtc_diag(pa_reac_data$sa5$mod$results)

#Nodesplit
pa_reac_data$sa5$nodesplit = mtc.nodesplit(pa_reac_data$sa5$mod$network)

summary(pa_reac_data$sa5$nodesplit)
#ANOHE
pa_reac_data$sa5$anohe = mtc.anohe(pa_reac_data$sa5$mod$network)

pdf("./figs/pain scales reactivity/pa_reac_anohe_sa5.pdf", height = 11, width = 8.5)
plot(summary(pa_reac_data$sa5$anohe))
dev.off()

summary(pa_reac_data$sa5$anohe$result.cons)


windows()
forest(relative.effect.table(pa_reac_data$sa5$mod$results),"drops")


# Add trials with imputed means to excluded list
imp_n = pa_reac_data$pa$gemtc$input %>% left_join(pa_reac[c("studlab","imputed_mean")] %>% distinct(), by = c("studlab"))  %>% filter(imputed_mean == "yes") %>% 
  summarise(sample = sum(n),
  treat = paste(treatment,collapse=' vs ')) %>%
  mutate(reason = "imputed mean")

pa_reac_excluded = rbind(pa_reac_excluded,imp_n)




#==========================================================================
# Power analysis
# Assumptions:
# 1 = I2 is 50%
# 2 = I2 is 70%
# Effect to detect = 2 points (1 MID)


#==========================================================================

#
pa_reac_sa1_power = power_table(data = pa_reac_data$sa1$mod$results)

chars_sa5 = netmeta_xl_chars(data = pa_reac %>% filter(!studlab %in% pa_reac_excluded$studlab),outcome = "pa_reac",ref = "drops",treat = "trt_group")

pa_reac_sa5_power = power_table(data = pa_reac_data$sa5$mod$results, direct_comps = chars_sa5$direct_zeros)

# Forest plot vs ref======


comp_names_sa1 = c("Sweet taste multisensory + TA","NNS + TA","Sweet taste + TA", "Placebo",
          "EBM multisensory + TA","Acetaminophen 30min + TA","Acetaminophen 60min + TA",
          "Sweet taste + N2O + TA","WFDRI + TA","Sweet taste alone","Repeated sweet taste",
          "Sweet taste + singing")


comp_names_sa5 = c("Sweet taste multisensory + TA","NNS + TA","Sweet taste + TA", "Placebo",
                   "EBM multisensory + TA","Acetaminophen 30min + TA",
                   "Sweet taste + N2O + TA","WFDRI + TA","Sweet taste alone","Repeated sweet taste",
                   "Sweet taste + singing")




momlinc_fp_p(data = pa_reac_sa5_power, names = comp_names_sa5, width = 7.5)





#================================== =
#================================== = 
#====== Sensitivity analysis plot====
#================================== =  
#================================== =

analyses = c("pa","sa1","sa2","sa3","sa4","sa5")
names = c("Primary \n analysis","Remove \n Boyle","PMA \n MR","RoB \n MR","Remove \n posters","Remove \n imputed mean")
treat_names = c("Sweet taste multisensory + TA","Sweet taste + TA","EBM multisensory + TA","Sweet taste + N2O + TA","Acetaminophen 60min + TA",
                "NNS + TA","Sweet taste alone","Repeated sweet taste","WFDRI + TA","Sweet taste + singing","Topical Anesthetic (TA)","Acetaminophen 30min + TA",
                "No treatment")


sa_table$Treatment = treat_names

reac_heatplot = rop_heatplots(data = pa_reac_data,analyses_index = analyses,pub_names = names, treats = treat_names)

windows()
reac_heatplot$plot


#============================================== =
####League table
#============================================== =
pa_reac_sa5names = c("Sweet taste \n multisensory + \n TA",
                     "Sweet taste + \n TA",
                     "Sweet taste + \n N2O + TA",
                     "EBM \n multisensory + \n TA",
                     "NNS + TA",
                     "Sweet taste \n alone",
                     "Repeated \n sweet taste",
                     "WFDRI + TA",
                     "Acetaminophen \n 30min + TA",
                     "Sweet taste + \n singing",
                     "Topical \n Anesthetic (TA)",
                     "No treatment")
reac_basicp = relative.effect(pa_reac_data$sa5$mod$results,t1 = c("drops"),preserve.extra = FALSE)
reac_results = as.data.frame(as.matrix(reac_basicp$samples)) %>% mutate(d.drops.drops = 0)
reac_order = pa_reac_data$sa5$mod$suc %>% mutate(pub_names = pa_reac_sa5names)

order = reac_order

windows()
league_plot(results = as.data.frame(as.matrix(reac_basicp$samples)) %>% mutate(d.drops.drops = 0), order = reac_order, textsize = 3.5)

#================================================
#Probability of 2 point or greater difference====
#================================================
prob_grt2 = relative.effect(pa_reac_data$sa5$mod$results,t1 = "drops_sweet", t2 = "drops_sweet_mult",preserve.extra = F)
summary(prob_grt2)
extracted = as.matrix(prob_grt2[[1]])


sum(extracted <= -2)/length(extracted)*100 + sum(extracted >= +2)/length(extracted)*100

#============================================== =
#Prob treatment mean < 6 points on the PIPP======
# Decision rule: Largest, lowest risk of bias during the procedure

reac_ma = pa_reac %>% filter(trt_group == "drops") %>% select(studlab,mean,sample_size,std_dev,actual_timepoint) %>% arrange(-sample_size) %>% mutate(std.err = std_dev/sqrt(sample_size)) %>%
  left_join(rop_data_study[c("studlab","oa_rob_sub")]) %>% filter(oa_rob_sub == "low") %>% arrange(std.err)

#============================================== =


baseline_pipp = reac_ma[1,]

#----Use mean and std err from above to create distribution of 
#baseline scores with uncertainty
post = as.data.frame(as.matrix(pa_reac_data$sa5$mod$results$samples))

n.sims = length(post$sd.d)

pipp_sim = rnorm(n.sims,baseline_pipp$mean,baseline_pipp$std.err)




#----Get basic parameters for outcomes that have two steps to get to drops (weird gemtc output quirk)
colnames(post)

reac_bp = post %>% select(-starts_with("d.drops."),-sd.d)
post = post %>% select(starts_with("d.drops."))

colnames(reac_bp)
basic_par = as.data.frame(as.matrix((relative.effect(pa_reac_data$sa5$mod$results, t1 = "drops",t2 = c("drops_phys","drops_N2O_sweet","sweet_rep","sweet_sing"),preserve.extra = F))$samples))

post = bind_cols(post,basic_par) 

#Baseline pipp + treatment effect
post = map_df(post,~pipp_sim + .) %>% add_column(drops = pipp_sim)

#Get probability less than 6 + quantiles
prob_nopain = purrr::map(post, ~(sum(. < 6)/n.sims)*100)

quantile(post$d.drops.drops_N2O_sweet, probs = c(0.025,0.5,0.975))


#Get babies with scores less than six
(absolute_graph = purrr::map(post,~quantile(.,c(0.025,0.5,0.975))) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    rownames_to_column("treatment") %>% gather(quantile, score,`2.5%`:`97.5%`) %>%
    mutate(below_six = round(pnorm(5.4,score,baseline_pipp$std_dev),2),
           six_thirteen = round(pnorm(c(12.4),score,baseline_pipp$std_dev) - below_six,2),
           high = 1-below_six - six_thirteen))

pipp_abs_graphs = function(data = absolute_graph, colours = c("#00CD00", "#FFD700", "#CD2626")){
  
  fifty = filter(data, quantile == "50%")
  lowcri = filter(data, quantile == "2.5%")
  highcri = filter(data, quantile == "97.5%")
  
  
  graphs = NULL
  data_graph = data %>% filter(quantile == "50%")
  
  for(i in seq_along(fifty[,1])){
    
    
    #Create Background scale
    t = data.frame(id = c(1,1,1),
                   pain = c("low","moderate","high"),
                   score = c(6,6,9))
    
    
    therm = t %>% ggplot(aes(x = id,y = score,fill = pain)) + geom_bar(stat = "identity") + 
      scale_fill_manual(values = rev(colours)) + 
      scale_x_discrete(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), breaks = seq(1,21, by = 2)) +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 22)) + guides(fill = FALSE) +
      
      geom_hline(yintercept = fifty[i,3][[1]], size = 2) +
      geom_hline(yintercept = lowcri[i,3][[1]], size = 1, linetype = "dashed") + 
      geom_hline(yintercept = highcri[i,3][[1]], size = 1, linetype = "dashed") +
      theme(plot.margin = unit(c(1.5,0,1.5,0.5),"cm"))
    
    #Create icon array
    person_data = list(low = data_graph[i,4][[1]], moderate = data_graph[i,5][[1]], high = data_graph[i,6][[1]])
    personograph(rev(person_data), plot.width = 1, 
                 colors = list(low = colours[1], moderate = colours[2], high = colours[3]),
                 icon.style = 6,
                 draw.legend = T, dimensions = c(5,20))
    a <- grid.grab()
    
    grid.arrange(therm,a,ncol =2, widths = c(1/8,7/8), top = textGrob(paste(fifty[i,1]), gp=gpar(fontsize = 15)))
    
    temp = grid.grab()
    
    graphs[[paste(filter(data,quantile == "50%")[i,1])]] = temp
    
  }
  graphs
}


#Select top three treat by sucra
graphs_pub = pa_reac_data$sa5$mod$suc %>% top_n(3,sucra) 


png("absolute_plot_panel.png", width= 20, height = 13, units = "in", res = 300)
windows()
grid.arrange(test_graphs[["drops"]],test_graphs[[graphs_pub[1,1]]],test_graphs[[graphs_pub[2,1]]],test_graphs[[graphs_pub[3,1]]], left = "PIPP score")
dev.off()
# #========================================================================================
# 
# 
# Outcome: Pain Recovery-----
# 
# #========================================================================================


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Prep Data
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

pa_recov_data = NULL

pa_recov_data$pa$gemtc = prep_gem(pa_recov)

pa_recov_data$pa$gemtc$data = pa_recov_data$pa$gemtc$data %>% left_join(rop_data_study[c("studlab","avg_pma","pub_type","oa_rob_sub")],by = c("study" = "studlab"))%>% 
  left_join(pa_recov[c("studlab","actual_timepoint")] %>% distinct(), by = c("study" = "studlab")) %>%
  mutate(actual_timepoint = ifelse(actual_timepoint == "1 min post-exam",1,0),
         oa_rob_sub = ifelse(oa_rob_sub == "low",0,1)) %>% replace_na(list(actual_timepoint = 0))


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Primary Analysis
# --------- Include crossovers
# --------- Mean difference outcome
# --------- Include imputed means
# --------- include scaled scores
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


pa_recov_data$pa$mod = set_net(pa_recov_data$pa$gemtc$data)

calc_n(data = pa_recov_data$pa$gemtc$data,input = pa_recov_data$pa$gemtc$input)


gemtc_diag(pa_recov_data$pa$mod$results)

summary(pa_recov_data$pa$mod$results)


#Analysis of heterogeneity and rough plot of effects versus reference 
# Nodeplit (primary)
pa_recov_nodesplit = mtc.nodesplit(pa_recov_data$pa$mod$network)
summary(pa_recov_nodesplit)

#ANOHE - Fits a UME, USE, and Cons model... good visualization
pa_recov_data$pa$anohe = mtc.anohe(pa_recov_data$pa$mod$network)

pdf("./figs/pain scales recovery/pa_recov_anohe.pdf", height = 11, width = 8.5)
plot(summary(pa_recov_data$pa$anohe))
gemtc::forest(relative.effect.table(pa_recov_data$pa$mod$results),"drops")
dev.off()


#Write league table to file
write.csv(pa_recov_data$pa$mod$league,"./tables/pa_recov_league.csv")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 1
# ----- Meta-regression on pma
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sa1_regressor = list(coefficient = "shared",
                     variable = "avg_pma",
                     control = "drops")

pa_recov_data$sa1$mod = set_net(data = pa_recov_data$pa$gemtc$data %>% drop_na(avg_pma),
                               model = "regression",regressor = sa1_regressor)



calc_n(data = pa_recov_data$pa$gemtc$data %>% drop_na(avg_pma),input = pa_recov_data$pa$gemtc$input)


gemtc_diag(pa_recov_data$sa1$mod$results)

summary(pa_recov_data$sa1$mod$results)

forest(relative.effect.table(pa_recov_data$sa1$mod$results),"drops")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 2
# ----- MR Risk of bias
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sa2_regressor = list(coefficient = "shared",
                     variable = "oa_rob_sub",
                     control = "drops")

pa_recov_data$sa2$mod = set_net(data = pa_recov_data$pa$gemtc$data,
                               model = "regression",regressor = sa2_regressor)

gemtc_diag(pa_recov_data$sa2$mod$results)

summary(pa_recov_data$sa2$mod$results)

forest(relative.effect.table(pa_recov_data$sa2$mod$results),"drops")


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
# Sensitivity 3
# ----- MR Actual timepoint (1 min = immediate)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

sa3_regressor = list(coefficient = "shared",
                     variable = "actual_timepoint",
                     control = "drops")

pa_recov_data$sa3$mod = set_net(data = pa_recov_data$pa$gemtc$data,
                                model = "regression",regressor = sa3_regressor)

gemtc_diag(pa_recov_data$sa3$mod$results)

summary(pa_recov_data$sa3$mod$results)

forest(relative.effect.table(pa_recov_data$sa3$mod$results),"drops")


#==========================================================================
# Power analysis
# Assumptions:
# Effect to detect = 2 points (1 MID)


#==========================================================================

#


pa_recov_power = power_table(data = pa_recov_data$pa$mod$results, direct_comps = recov_chars$direct_zeros)

# Forest plot vs ref======

recov_comp_names = c("NNS + TA","Sweet multisensory + TA","Sweet + TA",  "EBM multisensory + TA","Acetaminophen 30min + TA",
                     "Acetaminophen 60min + TA", "Morphine + TA", "NNS alone","No treatment", "Sweet taste alone")
                  


momlinc_fp_p(data = pa_recov_power, names = recov_comp_names, width = 7.5)


#============================================== =
####League table
#============================================== =
pa_recov_data$pa$mod$suc
pa_recov_panames = c("EBM \n multisensory + \n TA",
                     "Sweet taste \n multisensory + \n TA",
                     "NNS + TA",
                     "Morphine + TA",
                     "Sweet taste + \n TA",
                     "Acetaminophen \n 60min + TA",
                     "NNS alone",
                     "Acetaminophen \n 30min + TA",
                     "Topical \n Anesthetic (TA)",
                     "Sweet taste \n alone",
                     "No treatment")
recov_basicp = relative.effect(pa_recov_data$pa$mod$results,t1 = c("drops"),preserve.extra = FALSE)
recov_results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0)
recov_order = pa_recov_data$pa$mod$suc %>% mutate(pub_names = pa_recov_panames)


windows()
league_plot(results = as.data.frame(as.matrix(recov_basicp$samples)) %>% mutate(d.drops.drops = 0), order = recov_order, textsize = 3.5)

#================================== =
#================================== = 
#====== Sensitivity analysis plot====
#================================== =  
#================================== =

recov_analyses = c("pa","sa1","sa2","sa3")
recov_names = c("Primary \n analysis","PMA \n MR","RoB \n MR","Timepoint \n MR")

recov_suc_list = NULL

for(i in 1:length(recov_analyses)){
  recov_suc_list[[i]] = pa_recov_data[[recov_analyses[[i]]]]$mod$suc
  recov_suc_list[[i]] = recov_suc_list[[i]] %>% separate(treat, into = c("a","treat"), sep = "d.drops.") %>% select(-a)
}
recov_sa_table = recov_suc_list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="treat"), .)

colnames(recov_sa_table) = c("Treatment",recov_names)

treat_names = c("EBM multisensory + TA","Sweet multisensory + TA","NNS + TA", "Morphine + TA","Sweet + TA","Acetaminophen 60min + TA",
                "NNS alone", "Acetaminophen 30min + TA","Topical Anesthetic (TA)", "Sweet taste alone", "No treatment")

recov_sa_table$Treatment = treat_names
#Function starts

recov_sa_table_long = recov_sa_table %>% gather(variable,sucra,`Primary \n analysis`:`Timepoint \n MR`) %>% mutate(value2 = ifelse(!is.na(sucra),paste(round(sucra*100,0),"%",sep = ""),
                                                                                                                              paste(round(sucra*100,0)))) %>% 
  mutate(variable = as_factor(variable))

recov_sa_table_long$Treatment = factor(recov_sa_table$Treatment,
                                 levels = rev(recov_sa_table$Treatment))
windows()
recov_sa_table_long %>% ggplot(aes(y = Treatment, x = variable)) + 
  geom_tile(aes(fill = sucra), colour = "white", size = 2) + 
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


#================================================
#Probability of 2 point or greater difference====
#================================================
prob_grt2 = relative.effect(pa_recov_data$pa$mod$results,t1 = "drops_sweet_mult", t2 = "drops_ebm_mult",preserve.extra = F)
summary(prob_grt2)
extracted = as.matrix(prob_grt2[[1]])


sum(extracted <= -2)/length(extracted)*100 + sum(extracted >= +2)/length(extracted)*100

#============================================== =
#Prob treatment mean < 6 points on the PIPP======
# Selection rule: Baseline (TA) PIPP score in the largest/smallest se, low risk of bias study in the immediate phase
recov_ma = pa_recov %>% filter(trt_group == "drops") %>% select(studlab,mean,sample_size,std_dev,actual_timepoint) %>% arrange(-sample_size) %>% mutate(std.err = std_dev/sqrt(sample_size)) %>%
  left_join(rop_data_study[c("studlab","oa_rob_sub")]) %>%
  filter(oa_rob_sub == "low") %>% arrange(std.err)

recov_baseline_pipp = recov_ma[1,]
#============================================== =

#----Use mean and std err from above to create distribution of 
#baseline scores with uncertainty
recov_post = as.data.frame(as.matrix(pa_recov_data$pa$mod$results$samples))

n.sims = length(recov_post$sd.d)

recov_pipp_sim = rnorm(n.sims,recov_baseline_pipp$mean,recov_baseline_pipp$std.err)


#----Get basic parameters for outcomes that have two steps to get to drops (weird gemtc output quirk)
colnames(recov_post)

recov_bp = recov_post %>% select(-starts_with("d.drops."),-sd.d)
recov_post = recov_post %>% select(starts_with("d.drops."))

colnames(recov_bp)
recov_basic_par = as.data.frame(as.matrix((relative.effect(pa_recov_data$pa$mod$results, t1 = "drops",t2 = c("drops_ebm_mult","drops_phys","drops_sweet_mult","phys"),preserve.extra = F))$samples))

recov_post = bind_cols(recov_post,recov_basic_par) 

#Baseline pipp + treatment effect
recov_post = map_df(recov_post,~recov_pipp_sim + .) %>% add_column(drops = recov_pipp_sim)

summary(recov_post)
#Get probability less than 6 + quantiles
prob_nopain = purrr::map(recov_post, ~(sum(. < 6)/n.sims)*100)


#Get babies with scores less than six
(recov_absolute_graph = purrr::map(recov_post,~quantile(.,c(0.025,0.5,0.975))) %>% as.data.frame() %>% t() %>% as.data.frame() %>%
    rownames_to_column("treatment") %>% gather(quantile, score,`2.5%`:`97.5%`) %>%
    mutate(below_six = round(pnorm(5.4,score,recov_baseline_pipp$std_dev),2),
           six_thirteen = round(pnorm(c(12.4),score,recov_sd_pipp) - below_six,2),
           high = 1-below_six - six_thirteen))


recov_icon_array = pipp_abs_graphs(data = recov_absolute_graph)
  
graphs_pub = pa_recov_data$pa$mod$suc %>% top_n(3,sucra) 

windows()
grid.arrange(recov_icon_array[["drops"]],recov_icon_array[[graphs_pub[1,1]]],recov_icon_array[[graphs_pub[2,1]]],recov_icon_array[[graphs_pub[3,1]]], left = "PIPP score")


library(metafor)
library(fields)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(personograph)
source("./analyses/final/rop_import.R")
source("./functions/rankheat-plot-function.r")



rob = rop_data_study %>% select(studlab,rob_sg,rob_ac,rob_bp,rob_bo_ob,rob_bo_sub,rob_io,rob_sr,rob_other)



data = rob %>% filter(!(studlab %in% c("Zeraati 2015a","Zeraati 2015b","Zeraati 2015c"))) %>% gather(domain, `Risk of bias`, rob_sg:rob_other) %>% 
  mutate(domain = factor(domain, levels = c("rob_sg","rob_ac","rob_bp","rob_bo_ob","rob_bo_sub","rob_io","rob_sr","rob_other"),
                         labels = c("Sequence generation",
                                    "Allocation concealment",
                                    "Blinding of personnel",
                                    "Blinding of outcome assessors (objective)",
                                    "Blinding of outcome assessors (subjective)",
                                    "Incomplete outcome reporting",
                                    "Selective reporting",
                                    "Other")))
  
  
  data %>% ggplot(aes(y = studlab, x = domain, fill = `Risk of bias`)) + 
  geom_point(size = 6, shape = 21,stroke = 1.2) + scale_fill_manual(values = c("firebrick2","chartreuse3","yellow")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                             panel.background = element_blank()) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 1,vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank()) + scale_x_discrete(position = "top") 

ggsave("./figs/rob_table.pdf", height = 10, width = 5)


stacked = data %>% group_by(domain) %>% summarize(n = sum(!is.na(`Risk of bias`)),Low = sum(`Risk of bias` == "low", na.rm = TRUE)/n*100,
                                        Unclear = sum(`Risk of bias` == "unclear", na.rm = TRUE)/n*100,
                                        High = sum(`Risk of bias` == "high", na.rm = TRUE)/n*100) %>% gather(risk,prop,Low:High) %>% mutate(risk = factor(risk, levels = c("High","Unclear","Low")))


stacked %>% ggplot(aes(x = domain, y = prop, fill = risk)) + 
  geom_bar(stat = "identity", colour = "Black", size = 0.2) + scale_fill_manual(values = c("firebrick2","yellow","chartreuse3"),
                                                                    name = element_blank(),
                                                                    labels = c("High     ",
                                                                               "Unclear                                                               ",
                                                                               "Low                                                                   "),
                                                                    guide = guide_legend(reverse = TRUE)) + scale_x_discrete(limits = rev(levels(stacked$domain)),
                                                                                                                             expand = c(0,0.7)) + scale_y_continuous(limits = c(0,100),
                                                                                                                                                                                        expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.6,"line"),
        legend.text = element_text(size = 8),
        axis.line.x = element_line(colour = "black", size = 0.3),
        legend.position = "bottom",
        legend.justification = c(0,0),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.3), 
        legend.background = element_rect(size = 0.2, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0.3,0.6,0.3,0.3),"cm")) + coord_flip()

ggsave("./figs/rob_table_stacked.pdf", height = 2.5, width = 8.5)



#===========================================
# Excluded trials with reasons and n =======
#===========================================

excluded_from_ma = list(hr_reac = hr_reac_excluded,
                        hr_recov = hr_recov_excluded,
                        spo2_reac = os_reac_excluded,
                        spo2_recov = os_recov_excluded,
                        ae_reac = ae_reac_excluded)


#=============================================
# Rank heatplot ==============================
#=============================================
# install.packages("fields")
# install.packages("RColorBrewer")
# install.packages("circlize")



netheat_extract = function(sucra_loc = ae_reac$bugs$pa$rr,outcome = "adverse_events"){
 temp = as.data.frame(sucra_loc) %>% rownames_to_column() %>% select(rowname, `Mean SUCRA`)
 colnames(temp) = c("treatment",outcome)
 temp[outcome] = round(as.numeric(as.character(temp[[outcome]])),0)
 temp
}

pa_recov_sucra$pa_recov_sucra = round(pa_recov_sucra$pa_recov_sucra*100,0)

netheat_data = data.frame(pa_reac = round(pa_reac_sucra*100,0)) %>% rownames_to_column("treatment") %>%
  left_join(pa_recov_sucra, by = "treatment") %>% 
  left_join(netheat_extract(sucra_loc = ae_reac$bugs$pa$rr , outcome = "adverse_events"), by = "treatment") %>%
  left_join(netheat_extract(sucra_loc = cry_reac$bugs$pa$rr, outcome = "cry_time"), by = "treatment") %>%            
  left_join(netheat_extract(sucra_loc = hr_recov$bugs$pa$rr, outcome = "hr_recov"), by = "treatment") %>%  
  left_join(netheat_extract(sucra_loc = os_reac$bugs$pa$rr, outcome = "spo2_reac"), by = "treatment") %>%
  mutate(hr_reac = NA,
         spo2_recov = NA) %>% select(treatment,pa_reac,pa_recov_sucra,adverse_events,cry_time,hr_reac,hr_recov,
                                     spo2_reac,spo2_recov)

vector.outcome.names = c("PIPP reactivity",
                         "PIPP recovery",
                         "Adverse events",
                         "Crying time",
                         "Heart rate reactivity",
                         "Heart rate recovery",
                         "SpO2 reactivity",
                         "SpO2 recovery")


png("./figs/final/rop_heat_plot.png",res = 300,width = 2300,height = 2000)
rankheatplot(data = netheat_data,title.name ="Pain from ROP eye exams - Rank-heat plot based on SUCRA",
             cex = 0.65,
             pos.outcome.label = c(0.8,-0.1),
             show.numbers = TRUE,
             vector.outcomes = "vector.outcome.names")
dev.off()
chars

#================================================
#Probability of 2 point or greater difference====
#================================================
prob_grt2 = relative.effect(pa_reac_data$sa1$results,t1 = "drops_sweet", t2 = "drops_sweet_mult",preserve.extra = F)
summary(prob_grt2)
extracted = as.matrix(prob_grt2[[1]])


sum(extracted <= -2)/length(extracted)*100


recov_prob_grt2 = relative.effect(pa_recov_data$sa1$results,t1 = "drops_sweet_mult", t2 = "drops_ebm_mult",preserve.extra = F)
summary(recov_prob_grt2)
recov_extracted = as.matrix(recov_prob_grt2[[1]])


sum(recov_extracted <= -2)/length(extracted)*100

#================================================
#Outcome tables=================================
#================================================


sa_tables = function(data = pa_reac_data,names = c("pa","sa1","sa2","sa7"),n_wb = 1,pa_mr = F){
tables = NULL

 for(i in 1:(length(data)-n_wb)){
   results = summary(data[[i]]$results)
   rel_effects = t(as.data.frame(round(relative.effect.table(data[[i]]$results),2))[1,-(1:2)])
   rel_effects = as.data.frame(rel_effects) %>% rownames_to_column("treatment") %>% rename(cri = drops)
   
   sd = as.data.frame(results$summaries$quantiles) %>% rownames_to_column("treatment") %>% filter(treatment == "sd.d" | treatment == "B") %>% 
     select(treatment,`2.5%`,`50%`,`97.5%`) %>% mutate(cri = paste(round(`50%`,2)," (",round(`2.5%`,2)," to ",round(`97.5%`,2),")",sep = "")) %>% select(treatment,cri)
   
   bind = bind_rows(rel_effects,sd)
   
   tables[[i]] = as.data.frame(bind)
   
   names(tables)[[i]] = names[[i]]
 }

if(pa_mr == F){
tables[[1]] = tables[[1]] %>% add_row(treatment = "B", cri = NA)}

 for(i in (length(data)-n_wb+1):length(data)){
rel_effects = as.data.frame(data[[i]]$bugs$comp[grep("drops$",data[[i]]$bugs$comp[,1]),c(1:3)]) %>%
  rename(trt = `Comparison (Trt A vs. Trt B)`,
         mean = `Mean Difference of Trt A vs. Trt B`,
         cri = `95% CrI of Mean Difference`) %>% separate(trt,c("treatment","d"),sep = " vs. ") %>%
  mutate(cri = paste(mean," (",cri,")",sep = ""))  %>% select(-d,-mean)

sd = data.frame(treatment = c("sd.d","B"),
  cri = c(paste(round(data[[i]]$bugs$sd[5],2)," (",round(data[[i]]$bugs$sd[3],2)," to ",
           round(data[[i]]$bugs$sd[7],2),")",sep = ""),
          
          paste(round(data[[i]]$bugs$B[5],2)," (",round(data[[i]]$bugs$B[3],2)," to ",
                round(data[[i]]$bugs$B[7],2),")",sep = "")))
  
bind = bind_rows(rel_effects,sd)

tables[[i]] = as.data.frame(bind)

names(tables)[[i]] = names[[i]]
}

#tables = tables %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="treatment",suffix = c(".1",".2")), .)

#colnames(tables) = c("treatment",names)



tables
}

t = sa_tables()

#================================================
#Prob treatment mean < 6 points on the PIPP======
#================================================

#----Meta-analysis of baseline mean
ma = pa_reac %>% filter(trt_group == "drops") %>% select(studlab,mean,sample_size,std_dev) %>% arrange(-sample_size) %>% mutate(std.err = std_dev/sqrt(sample_size))

#If using meta-analysis (hard to justify given heterogeneity)
# baseline_pipp = rma(yi = mean,vi = std.err, data= ma)

#If using largest study
baseline_pipp = ma %>% filter(sample_size == max(sample_size)) %>% rename(se = std.err,
                                                                          b = mean)

#----Use mean and std err from above to create distribution of 
#baseline scores with uncertainty
post = as.data.frame(as.matrix(pa_reac_data$pa$results$samples))

n.sims = length(post$sd.d)

pipp_sim = rnorm(n.sims,baseline_pipp$b[[1]],baseline_pipp$se)


#----Get basic parameters for outcomes that have two steps to get to drops (weird gemtc output quirk)
post = post %>% select(-d.drops_sweet.drops_N2O_sweet,-d.placebo.sweet_rep,-d.placebo.sweet_sing,-sd.d)

basic_par = as.data.frame(as.matrix((relative.effect(pa_reac_data$pa$results, t1 = "drops",t2 = c("drops_N2O_sweet","sweet_rep","sweet_sing"),preserve.extra = F))$samples))

post = bind_cols(post,basic_par)

#Baseline pipp + treatment effect
post = map_df(post,function(x) pipp_sim + x)

#Negative numbers not allowed
post[post < 0] = 0

#Get probability less than 6 + quantiles
prob_nopain = purrr::map(post, function(x) (sum(x < 6)/n.sims)*100)





#Get babies with scores less than six
absolute_graph = purrr::map(post,base::mean) %>% as.data.frame() %>% gather(treatment, score) %>%
  mutate(below_six = pnorm(6,score,baseline_pipp$std_dev),
         six_thirteen = pnorm(c(13),score,baseline_pipp$std_dev) - below_six,
         high = 1-below_six - six_thirteen)


data_graph = purrr::map(absolute_graph[,3:5],as.numeric) 
treats = absolute_graph[,1]

test = list(low = data_graph[[1]][1], moderate = data_graph[[2]][1], high = data_graph[[3]][1])
personograph(test, colors = list(low = colours[1], moderate = colours[2], high = colours[3]))

colours = c("#00CD00", "#FFD700", "#CD2626")


absolute_graph %>% slice(1) %>% select(treatment,score)

t = data.frame(id = c("pipp","pipp","pipp"),
           pain = c("low","moderate","high"),
           score = c(6,13,22))

t %>% ggplot(aes(y = score,colour = pain)) + geom_bar()

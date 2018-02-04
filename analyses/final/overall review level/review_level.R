# library(metafor)
library(fields)
library(RColorBrewer)
library(circlize)
library(tidyverse) #fields loads a map function which interferes with this
library(personograph)
library(gridExtra)
library(gridBase)
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


#===========================================
# Additional explainations           =======
#===========================================
explanations = tibble(trial = c("Kleberg 2008"),
                      note = c("GA,BW,PMA presented as medians by site. Mean of medians used"))


#=============================================
# Rank heatplot ==============================
#=============================================

#HR Reac needs to be row of NAs
#OS recov needs to be row of NAs



pa_reac_data$sa5$mod$suc = pa_reac_data$sa5$mod$suc %>% add_row(treat = c("d.drops.phys","d.drops.drops_acet60", "d.drops.drops_morph"),
                                     sucra = c(NA,NA,NA))

all_out_list = list(pa_reac = pa_reac_data$sa5$mod$suc, pa_recov = pa_recov_data$pa$mod$suc,
                    hr_recov = hr_recov_data$pa$mod$suc, os_reac = os_sucra,
                    cry = cry_reac_data$pa$mod$suc, ae_reac = ae_sucra, ae_recov = recov_ae_sucra)

names = c("Pain \n reactivity","Pain \n regulation",
          "Heart Rate \n regulation",
          "SpO2 \n reactivity",
          "Cry",
          "Adverse events \n reactivity",
          "Adverse events \n regulation")


order =  c("Treatment",
           "Pain \n reactivity","Pain \n regulation",
           "Heart Rate \n reactivity",
           "Heart Rate \n regulation",
           "SpO2 \n reactivity",
           "SpO2 \n regulation",
           "Cry",
           "Adverse events \n reactivity",
           "Adverse events \n regulation")

treat_names = c("Sweet taste multisensory + TA",
                "Sweet taste + TA",
                "Sweet taste + N2O + TA",
                "EBM multisensory + TA",
                "NNS + TA",
                "Sweet taste alone",
                "Repeated sweet taste",
                "WFDRI + TA",
                "Acetaminophen 30min + TA",
                "Sweet taste + singing",
                "Topical Anesthetic (TA)",
                "No treatment",
                "NNS alone",
                "Acetaminophen 60min + TA",
                "Morphine + TA")


all_out_heat = all_out_list %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="treat"), .)
  
names(all_out_heat) = c("Treatment",names)

all_out_heat = all_out_heat %>% mutate(`Heart Rate \n reactivity` = NA,
                        `SpO2 \n regulation` = NA)

all_out_heat$Treatment = treat_names

windows()
heatplot(all_out_heat %>% select(order))
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


##### Study characteristics

(n_reports = length(rop_data_study$studlab))
studies = n_reports - as.numeric(rop_data_study %>% filter(connected == "yes") %>% summarise(n()-1))

(rcts = (rop_data_study %>% filter(design == "Parallel") %>% select(studlab))[[1]] %>% length())
(xover = (rop_data_study %>% filter(design == "Crossover") %>% select(studlab))[[1]] %>% length())

treats = rop_data_arm %>% select(studlab,treatment,trt_group)
trtnames = treats$trt_group %>% unique()

proc = rop_data_study %>% select(studlab,speculum,scleral_dep,con_swad)

summary(proc)
proc %>% filter(con_swad == "Swaddle")

study_chars = rop_data_study %>% select(studlab,design,pub_type,control:trt3,method,speculum,scleral_dep,avg_pma,avg_bw) %>%
  mutate(trts = ifelse(!is.na(trt3),paste(control,"vs",trt1,"vs",trt2,"vs",trt3,sep = " "),
                       ifelse(!is.na(trt2),paste(control,"vs",trt1,"vs",trt2,sep = " "),
                              paste(control,"vs",trt1,sep = " ")
                              )
                       )
         ) %>% select(-c(trt1:trt3))
write.csv(file = "./tables/final/rop_studychars.csv",study_chars)



####====PRISMA Included in meta-analysis currently need to subtract 1 because
# os_reac done quickly and includes zeraati 2015 instead of 2016
length(c(as.character(pa_reac$studlab),
  as.character(pa_recov$studlab),
  
  as.character(ae_reac$data$studlab),
  as.character(cry_reac$data$studlab),
  
  as.character(os_reac$data$studlab),
  as.character(os_recov$data$studlab),
  
  as.character(hr_reac$data$studlab),
  as.character(hr_recov$data$studlab)) %>% unique()
)


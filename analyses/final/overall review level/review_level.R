source("./analyses/final/rop_import.R")

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
# Excluded trials with reasons and n
#===========================================

excluded_from_ma = list(hr_reac = hr_reac_excluded,
                        hr_recov = hr_recov_excluded,
                        spo2_reac = os_reac_excluded,
                        spo2_recov = os_recov_excluded,
                        ae_reac = ae_reac_excluded)


#=============================================
# Rank heatplot
#=============================================
# install.packages("fields")
# install.packages("RColorBrewer")
# install.packages("circlize")

library(fields)
library(RColorBrewer)
library(circlize)
source("./functions/rankheat-plot-function.r")

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


#Probability of 2 point or greater difference
prob_grt2 = relative.effect(pa_reac_data$sa1$results,t1 = "drops_sweet", t2 = "drops_sweet_mult",preserve.extra = F)
summary(prob_grt2)
extracted = as.matrix(prob_grt2[[1]])


sum(extracted <= -2)/length(extracted)*100


recov_prob_grt2 = relative.effect(pa_recov_data$sa1$results,t1 = "drops_sweet_mult", t2 = "drops_ebm_mult",preserve.extra = F)
summary(recov_prob_grt2)
recov_extracted = as.matrix(recov_prob_grt2[[1]])


sum(recov_extracted <= -2)/length(extracted)*100


#Outcome tables

results = summary(pa_reac_data$pa$results)
rel_effects = t(as.data.frame(round(relative.effect.table(pa_reac_data$pa$results),2))[1,-(1:2)] )
rel_effects = as.data.frame(rel_effects) %>% rownames_to_column("treatment") %>% rename(pa = drops)


as.data.frame(results$summaries$quantiles) %>% rownames_to_column("treatment") %>% filter(treatment == "sd.d") %>% 
  select(treatment,`2.5%`,`50%`,`97.5%`) %>% mutate(cri = paste(round(`50%`,2)," (",round(`2.5%`,2)," to ",round(`97.5%`,2),")",sep = "")) %>% select(treatment,cri)

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

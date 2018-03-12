library(grid)

pa_reac_sa5names_fp = c("Sweet taste multisensory + TA",
                     "Sweet taste + TA",
                     "Sweet taste + N2O + TA",
                     "EBM multisensory + TA",
                     "NNS + TA",
                     "Sweet taste alone",
                     "Acetaminophen 30min + TA",
                     "WFDRI + TA",
                     "Repeated sweet taste",
                     "Sweet taste + singing",
                     "No treatment")

df = as.data.frame(as.matrix(reac_basicp$samples)) %>% map_df(~quantile(.,probs = c(0.025,0.5,0.975))) %>% 
  gather(comp,value) %>% mutate(est = rep(c("lower","mean","upper"),11)) %>% 
  spread(est,value) %>% arrange(mean) %>% mutate(meandiff = paste(round(mean,2)," (",round(lower,2)," to ",round(upper,2),")",sep = ""),
                                               comp_n = pa_reac_sa5names_fp)



plt_data = df %>% select(mean,lower,upper) 

table = df %>% select(comp_n,meandiff)

windows()
forestplot(rbind(c("Comparison","Mean Difference (95% CrI)"),table),
           c(NA,plt_data$mean),
           c(NA,plt_data$lower),
           c(NA,plt_data$upper),
           graph.pos = 2,
           is.summary = c(TRUE,rep(FALSE,11)),
           hrzl_lines = gpar(col="#444444"),
           align = c("l",rep("c",5),"l"),
           boxsize = 0.4,
           colgap = unit(3,"mm"),
           lineheight = unit(7,"mm"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 11, family = "calibri"),
                            ticks = gpar(fontsize = 20, family = "calibri")),
           col = fpColors(box = "mediumpurple", line = "midnightblue"),
           xlab = "Compared to drops alone",
           title = "PIPP Reactivity (Intervention vs Anesthetic eye drops alone)"
           )


  
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
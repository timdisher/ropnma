

#======================
# Desired workflow
#======================

test = prep_wb(pa_reac)
test$ready = test$arm_wide %>% mutate(
  se_1 = sd_1/sqrt(n_1),
  se_2 = sd_2/sqrt(n_2),
  se_3 = sd_3/sqrt(n_3),
  se_4 = sd_4/sqrt(n_4)) %>% select(matches("t_"),matches("y_"),matches("se_"),na) %>% arrange(na)

(wb_test = nma_winbugs_datalist(test$ready,test$treatments,contrast = FALSE))

wb_test$mx = as.vector(test$ready %>% filter(t_1 == 1) %>% summarise(mx = mean(y_1)))[[1]]


metaregtest_armwise = bugs(wb_test,NULL,params_mr,model.file = MODELFILE.re_armdata_meta,
                   n.chains = 3, n.iter = 100000, n.burnin = 40000, n.thin = 10,
                   bugs.directory = bugsdir, debug = F)

loo = as.data.frame(metaregtest_armwise$summary)


test$ready %>% filter(t_1 == 1) %>% mutate(y2_diff = y_2 - y_1,
                                           y3_diff = y_3 - y_1,
                                           y4_diff = y_4 - y_1) %>% gather(diff,value,y2_diff:y4_diff) %>% gather(trt,num, t_2:t_4) %>%
  select(y_1,value,num) %>% na.omit %>% left_join(test$treatments, by = c("num" = "t")) %>%
  ggplot(aes(y = value, x = y_1)) + geom_point() + geom_smooth(method = "lm", se = FALSE,colour = "black") + facet_wrap(~description)



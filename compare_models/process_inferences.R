require(tidyverse)
require(cowplot)
require(purrr)
require(broom)

theme_set(theme_classic() + theme(axis.title = element_text(size=11), axis.text = element_text(size=10)))

path = "inference/"

################################# 
########### model fit ###########

model.fit <- read.csv(paste0(path,"model_fits.csv"))
model.fit %>% 
    group_by(dataset) %>%
    arrange(dataset, AICc) %>%
    mutate(AICc_rank = rank(AICc)) %>% 
    mutate(deltaAICc = min(AICc) - AICc) -> model.aicc
    
model.aicc %>% 
    filter(AICc_rank == 1) %>% 
    group_by(model) %>%
    tally() 
# 1  JTT-G     4
# 2   LG-G     1
# 3  WAG-G    85
# 4 WAG-No    10


    
model.aicc %>% 
    filter(AICc_rank == 8) %>%
    ungroup() %>%
    group_by(model) %>%
    tally() 
# 1  JC69-G    33
# 2 JC69-No    44
   
model.aicc %>% 
    filter(AICc_rank == 7) %>%
    ungroup() %>%
    group_by(model) %>%
    tally()      
# 1  JC69-G    44
# 2 JC69-No    39
# 3   JTT-G     6
# 4  JTT-No     1
# 5    LG-G     3
# 6   LG-No     2


model.aicc %>% 
    ungroup() %>%
    filter(model %in% c("JC69-No", "JC69-G")) %>%
    group_by(AICc_rank) %>%
    tally()  
#   AICc_rank     n
#       <dbl> <int>
# 1         4    10
# 2         5     2
# 3         6     8
# 4         7    83
# 5         8    77  

model.aicc %>% 
    filter(AICc_rank == 2) %>% 
    ungroup() %>%
    summarize(mean(deltaAICc), min(abs(deltaAICc)), max(abs(deltaAICc)))
# # A tibble: 1 x 3
#   `mean(deltaAICc)` `min(abs(deltaAICc))` `max(abs(deltaAICc))`
#                    <dbl>                 <dbl>                 <dbl>
# 1                -1958.132              27.5292               29917.6



model.aicc %>% filter(AICc_rank <= 6) %>%
    ggplot(aes(x = AICc_rank, y = deltaAICc, group = dataset)) + geom_line()





path = "inference/"
### Read in hyphy inferences ###
files <- dir(path = path, pattern = "*csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(path, .)))) %>%
  unnest() %>%
  separate(filename, c("dataset", "blahh", "model", "blahhhh", "blahhhhhhh"), sep="\\.") %>%
  dplyr::select(dataset, model, site, rate) -> rawrates


#### No normalization
rates <- rawrates %>% 
            group_by(dataset, model) %>% 
            spread(model, rate)

### Normalized by mean   
normmean <- ratesraw %>% 
                group_by(dataset, model) %>% 
                mutate(normrate = rate/mean(rate)) %>%
                select(-rate) %>%
                spread(model, normrate)

## Normalized by median
normmedian <- ratesraw %>% 
                group_by(dataset, model) %>% 
                mutate(normrate = rate/median(rate)) %>%
                select(-rate) %>%
                spread(model, normrate) 

## z score-ify
zrates <- ratesraw %>% 
                group_by(dataset, model) %>% 
                mutate(zscore = (rate - mean(rate))/sd(rate)) %>%
                select(-rate) %>%
                spread(model, zscore) 

# 
# ### Normalized by mean of non-zero/non-infinite
# ratesraw %>% 
#     group_by(dataset, model) %>% 
#     filter(rate >1e-8, rate < 100) %>%
#     summarize(mean_noextremes = mean(rate), med_noextremes = median(rate)) %>% 
#     right_join(ratesraw) -> raw.withnoextreme
#         
#      
# ### Normalized by non-extreme mean   
# normmean.noextreme <- raw.withnoextreme %>% 
#                         group_by(dataset, model) %>% 
#                         mutate(normrate = rate/mean_noextremes) %>%
#                         select(-rate,-mean_noextremes, -med_noextremes) %>%
#                         spread(model, normrate)
# 
# ## Normalized by non-extreme median   
# normmedian.noextreme <- raw.withnoextreme %>% 
#                         group_by(dataset, model) %>% 
#                         mutate(normrate = rate/med_noextremes) %>%
#                         select(-rate) %>%
#                         spread(model, normrate) 
#                 

      
###########################################################################################       
## Optimizing branch lengths with a gamma has little to no influence on rates, with exception of hard-to-estimate sites


zrates %>%
    group_by(dataset) %>%
    do(corLG = cor(.$`LG-G`, .$`LG-No`, method = "spearman"),
       corWAG = cor(.$`WAG-G`, .$`WAG-No`, method = "spearman"),
       corJTT = cor(.$`JTT-G`, .$`JTT-No`, method = "spearman"),
       corJC69 = cor(.$`JC69-G`, .$`JC69-No`, method = "spearman")) %>%
    unnest(corLG, corWAG, corJTT, corJC69) %>%
    summarize(mean(corLG), 
              mean(corWAG), 
              mean(corJTT), 
              mean(corJC69), 
              sd(corLG), 
              sd(corWAG), 
              sd(corJTT), 
              sd(corJC69)) %>% print.data.frame()
#   mean(corLG) mean(corWAG) mean(corJTT) mean(corJC69)   sd(corLG)  sd(corWAG)
# 1   0.9983799    0.9983811    0.9984812     0.9978295 0.001992638 0.002258409
#   sd(corJTT) sd(corJC69)
# 1 0.00183936 0.002369673


rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    do(corLG = cor(.$`LG-G`, .$`LG-No`, method = "spearman"),
       corWAG = cor(.$`WAG-G`, .$`WAG-No`, method = "spearman"),
       corJTT = cor(.$`JTT-G`, .$`JTT-No`, method = "spearman"),
       corJC69 = cor(.$`JC69-G`, .$`JC69-No`, method = "spearman")) %>%
    unnest(corLG, corWAG, corJTT, corJC69) %>%
    summarize(mean(corLG), 
              mean(corWAG), 
              mean(corJTT), 
              mean(corJC69), 
              sd(corLG), 
              sd(corWAG), 
              sd(corJTT), 
              sd(corJC69)) %>% print.data.frame()
#   mean(corLG) mean(corWAG) mean(corJTT) mean(corJC69)   sd(corLG)  sd(corWAG)
# 1   0.9986282    0.9986402    0.9986912     0.9977772 0.001514524 0.001696184
#    sd(corJTT) sd(corJC69)
# 1 0.001386472 0.002517279

rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    mutate(rankLGG = rank(`LG-G`), rankLGNo = rank(`LG-No`)) %>%
    do(fitLG = lm(rankLGG ~ rankLGNo, data=.)) %>%
    tidy(fitLG) %>%
    filter(term == "(Intercept)", p.value <= 0.01) %>%
    nrow()   
[1] 0    



rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    mutate(rankWAGG = rank(`WAG-G`), rankWAGNo = rank(`WAG-No`)) %>%
    do(fitWAG = lm(rankWAGG ~ rankWAGNo, data=.)) %>%
    tidy(fitWAG) %>%
    filter(term == "(Intercept)", p.value <= 0.01) %>%
    nrow()   
[1] 0    



rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    mutate(rankJTTG = rank(`JTT-G`), rankJTTNo = rank(`JTT-No`)) %>%
    do(fitJTT = lm(rankJTTG ~ rankJTTNo, data=.)) %>%
    tidy(fitJTT) %>%
    filter(term == "(Intercept)", p.value <= 0.01) %>%
    nrow()   
[1] 0    


rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    mutate(rankJC69G = rank(`JC69-G`), rankJC69No = rank(`JC69-No`)) %>%
    do(fitJC69 = lm(rankJC69G ~ rankJC69No, data=.)) %>%
    tidy(fitJC69) %>%
    filter(term == "(Intercept)", p.value <= 0.01) %>%
    nrow()   
[1] 0    









####################################################################################################


##### NOW compare across models ####

rates %>%
    group_by(dataset) %>%
    na.omit() %>%
    do(LGWAG = cor(.$`LG-No`, .$`WAG-No`, method = "spearman"),
       LGJTT = cor(.$`LG-G`, .$`JTT-No`, method = "spearman"),
       LGJC = cor(.$`LG-No`, .$`JC69-No`, method = "spearman"),
       WAGJTT = cor(.$`WAG-No`, .$`JTT-No`, method = "spearman"),
       WAGJC = cor(.$`LG-No`, .$`JC69-No`, method = "spearman"),
       JTTJC = cor(.$`JTT-No`, .$`JC69-No`, method = "spearman")) %>%
    unnest(LGWAG, LGJTT, LGJC, WAGJTT, WAGJC, JTTJC) %>% 
    summarize(mean(LGWAG), 
              mean(LGJTT), 
              mean(LGJC), 
              mean(WAGJTT), 
              mean(WAGJC), 
              mean(JTTJC), 
              sd(LGWAG), 
              sd(LGJTT), 
              sd(LGJC), 
              sd(WAGJTT), 
              sd(WAGJC), 
              sd(JTTJC)) %>% print.data.frame()
#   mean(LGWAG) mean(LGJTT) mean(LGJC) mean(WAGJTT) mean(WAGJC) mean(JTTJC)
# 1    0.992614   0.9877326     0.9507    0.9946388      0.9507   0.9571543
#     sd(LGWAG)   sd(LGJTT)   sd(LGJC) sd(WAGJTT)  sd(WAGJC)  sd(JTTJC)
# 1 0.003328295 0.005512693 0.01947251 0.00260622 0.01947251 0.01906305




#####################################################################################










       
       
       
       
       
       
       
       
       
       
       
       
       # 
#        
#                 
#                 
#                 ratesraw %>% filter(dataset != "1DD8_A") %>% spread(model, rate) %>% gather(model, rate, `JC69-Gamma`:`WAG-Gamma`) %>% group_by(model,dataset) %>% do(fit = lm(rate ~ `WAG-No`, data=.)) -> dffits
#                 
# ggplot(rates.nonorm, aes(x = `JTT-Gamma`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ facet_grid(~dataset) -> a
# ggplot(rates.nonorm, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()+ facet_grid(~dataset) -> b
#                 
# ggplot(normmedian.noextreme, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()+ facet_grid(~dataset) -> d
# ggplot(normmean.noextreme, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()+ facet_grid(~dataset) -> e
# 
# 
# 
# plot_grid(a,b,d,e,nrow=4)
#           
#                 
#                 
#                 
#                 
#                 
# 
# 
# 
#                 
# ### We can also normalize by mean and median of non-zero rates
# ### We should probably do z-scores also
# 
# ggplot(rates.nonorm, aes(x = `LG-Gamma`, y = `LG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG G vs noRV")+ scale_x_log10() + scale_y_log10() + facet_grid(~dataset)
# ggplot(rates.nonorm, aes(x = `WAG-Gamma`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG G vs noRV")+ scale_x_log10() + scale_y_log10() + facet_grid(~dataset)
# ggplot(rates.nonorm, aes(x = `JTT-Gamma`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG G vs noRV")+ scale_x_log10() + scale_y_log10() + facet_grid(~dataset)
# ggplot(rates.nonorm, aes(x = `JC69-Gamma`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JC69 G vs noRV") + scale_x_log10() + scale_y_log10()+ facet_grid(~dataset)
# 
# 
# ggplot(zrates, aes(x = `LG-No`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-WAG") + xlim(-2,2) + ylim(-2,2) + facet_grid(~dataset) -> a
# ggplot(zrates, aes(x = `JTT-No`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JTT-WAG")+ xlim(-2,2) + ylim(-2,2) +facet_grid(~dataset) -> b
# ggplot(zrates, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ xlim(-2,2) + ylim(-2,2) + facet_grid(~dataset) -> d
# ggplot(zrates, aes(x = `LG-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JC69")+ xlim(-2,2) + ylim(-2,2) +facet_grid(~dataset) -> e
# plot_grid(a,b,d,e,nrow=4)
# 
# 
# ggplot(normmedian, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()+ facet_grid(~dataset) -> a
# ggplot(normmean, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()+ facet_grid(~dataset) -> b
# plot_grid(a,b,nrow=2)
# 
# 
# ggplot(aa.spread, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()
# ggplot(aa.spread, aes(x = `LG-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JC69")+ scale_x_log10() + scale_y_log10()
# ggplot(aa.spread, aes(x = `WAG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("WAG-JTT")+ scale_x_log10() + scale_y_log10()
# ggplot(aa.spread, aes(x = `WAG-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("WAG-JC69")+ scale_x_log10() + scale_y_log10()
# ggplot(aa.spread, aes(x = `JTT-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JTT-JC69")+ scale_x_log10() + scale_y_log10()
# 
# 
# 
# ##################### AA plots ##########################
# paa1 <- ggplot(aa.spread, aes(x = `LG-Gamma`, y = `LG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG G vs noRV")+ scale_x_log10() + scale_y_log10() + facet_grid(~dataset)
# paa2 <- ggplot(aa.spread, aes(x = `WAG-Gamma`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("WAG G vs noRV") + scale_x_log10() + scale_y_log10()
# paa3 <- ggplot(aa.spread, aes(x = `JTT-Gamma`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JTT G vs noRV") + scale_x_log10() + scale_y_log10()
# paa4 <- ggplot(aa.spread, aes(x = `JC69-Gamma`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JC69 G vs noRV") + scale_x_log10() + scale_y_log10()
# plot_grid(paa1, paa2, paa3, paa4, nrow=1)
# 
# paa4 <- ggplot(aa.spread, aes(x = `LG-No`, y = `WAG-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-WAG")+ scale_x_log10() + scale_y_log10()
# paa5 <- ggplot(aa.spread, aes(x = `LG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JTT")+ scale_x_log10() + scale_y_log10()
# paa6 <- ggplot(aa.spread, aes(x = `LG-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("LG-JC69")+ scale_x_log10() + scale_y_log10()
# paa7 <- ggplot(aa.spread, aes(x = `WAG-No`, y = `JTT-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("WAG-JTT")+ scale_x_log10() + scale_y_log10()
# paa8 <- ggplot(aa.spread, aes(x = `WAG-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("WAG-JC69")+ scale_x_log10() + scale_y_log10()
# paa9 <- ggplot(aa.spread, aes(x = `JTT-No`, y = `JC69-No`)) + geom_point() + geom_abline(color="blue") + ggtitle("JTT-JC69")+ scale_x_log10() + scale_y_log10()
# 
# plot_grid( paa4, paa5, paa6, paa7, paa8, paa9, nrow=1)
# 
# plot_grid(paa1, paa2, paa3,, nrow=2) -> aa.plots
# 
# 
# 
# #################### NUC plots #############################
# 
# pnuc1 <- ggplot(nuc.spread1, aes(x = GTR_G, y = GTR_NoRV)) + geom_point() + geom_abline(color="blue") + ggtitle("GTR G vs noRV")
# pnuc2 <- ggplot(nuc.spread1, aes(x = HKY85_G, y = HKY85_NoRV)) + geom_point() + geom_abline(color="blue") + ggtitle("HKY85 G vs noRV")
# pnuc3 <- ggplot(nuc.spread1, aes(x = JC69_G, y = JC69_NoRV)) + geom_point() + geom_abline(color="blue") + ggtitle("JC69 G vs noRV")
# 
# pnuc4 <- ggplot(nuc.spread1, aes(x = GTR_G, y = HKY85_G)) + geom_point() + geom_abline(color="blue") + ggtitle("GTR-HKY")
# pnuc5 <- ggplot(nuc.spread1, aes(x = GTR_G, y = JC69_G)) + geom_point() + geom_abline(color="blue") + ggtitle("GTR-JC69")
# pnuc6 <- ggplot(nuc.spread1, aes(x = HKY85_G, y = JC69_G)) + geom_point() + geom_abline(color="blue") + ggtitle("HKY-JC69")
# 
# 
# plot_grid(pnuc1, pnuc2, pnuc3, pnuc4, pnuc5, pnuc6, nrow=2) -> nuc.plots
# 





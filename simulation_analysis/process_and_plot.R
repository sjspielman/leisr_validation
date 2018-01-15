require(tidyverse)
require(cowplot)
require(broom)

theme_set(theme_classic() + theme(axis.title = element_text(size=11), axis.text = element_text(size=10)))
inpath = "simulation_rates/"
figpath = "../figures/"


#############################################################################################
#################################### Read and process data###################################
#############################################################################################

### Read in hyphy inferences ###
files <- dir(path = inpath, pattern = "*hyphy.csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(inpath, .)))) %>%
  unnest() %>%
  separate(filename, c("sim", "repl", "alg", "hyphy"), sep="_") %>%
  unite("method", c(alg, hyphy)) %>%
  separate(method, c("method", "byebye"), sep = "\\.") %>% 
  dplyr::select(sim, repl, site, rate, method)  -> hyphy.rates
  
### Read in r4s inferences ###
files <- dir(path = inpath, pattern = "*r4s.csv")
data_frame(filename = files) %>%
  mutate(file_contents = map(filename,
           ~ read_csv(file.path(inpath, .)))) %>%
  unnest() %>%
  separate(filename, c("sim", "repl", "alg", "het", "blah"), sep="_") %>%
  unite("method",c("alg", "het")) %>%
  dplyr::select(sim, repl, method, site, rate) -> r4s.rates

### Merge to a single data frame of rates ###
full.rates <- rbind(hyphy.rates, r4s.rates)
full.rates %>% 
    group_by(method, sim, repl) %>%
    mutate(normrate = rate / mean(rate)) -> full.rates
full.rates %>%
    dplyr::select(-rate) %>%
    spread(method, normrate) -> full.spread

#############################################################################################
###################################### Generate barplot #####################################
#############################################################################################

full.spread %>%
    group_by(sim, repl) %>%
    do(tidy(cor.test(.$homo_hyphy, .$ML_homo))) %>%
    mutate(r2 = estimate**2, comp = "homo_homo") %>%
    select(sim, repl, r2, comp) -> r1

full.spread %>%
    group_by(sim, repl) %>%
    do(tidy(cor.test(.$homo_hyphy, .$ML_gamma))) %>%
    mutate(r2 = estimate**2, comp = "homo_gamma") %>%
    select(sim, repl, r2, comp) -> r2

full.spread %>%
    group_by(sim, repl) %>%
    do(tidy(cor.test(.$gamma_hyphy, .$ML_homo))) %>%
    mutate(r2 = estimate**2, comp = "gamma_homo") %>%
    select(sim, repl, r2, comp)  -> r3
    
full.spread %>%
    group_by(sim, repl) %>%
    do(tidy(cor.test(.$gamma_hyphy, .$ML_gamma))) %>%
    mutate(r2 = estimate**2, comp = "gamma_gamma") %>%
    select(sim, repl, r2, comp)  -> r4


rsquareds <- rbind(r1,r2) %>% rbind(r3) %>% rbind(r4)

rsquareds %>% 
    group_by(sim, comp) %>%
    summarize(meanr2 = mean(r2), se = sd(r2)/sqrt(n())) %>%
    mutate(lower = meanr2 -se/2, upper = meanr2 + se/2) -> rsquareds.summary


rsquareds.summary$sim <- factor(rsquareds.summary$sim, levels=c("sim25taxa", "sim50taxa", "sim100taxa"), labels=c("25 taxa", "50 taxa", "100 taxa"))
rsquareds.summary$comp <- factor(rsquareds.summary$comp, levels=c("homo_homo", "gamma_gamma", "homo_gamma", "gamma_homo"), labels=c("LEISR ~ R4S", "LEISR+G ~ R4S+G", "LEISR ~ R4S+G", "LEISR+G ~ R4S"))    


dodge <- position_dodge(width=0.9)
r2.bars <- ggplot(rsquareds.summary, aes(x = sim, y = meanr2, fill = comp)) +                  
                  geom_bar(stat="identity", position=dodge) + 
                  geom_errorbar(aes(ymin=lower, ymax=upper), position = dodge, width=0.25) + 
                  ylab(expression(R^{2})) + xlab("Simulation Set") + 
                  coord_cartesian(ylim = c(0.85, 1)) +
                  scale_fill_brewer(palette = "Set1", name = "Comparison")
save_plot(paste0(figpath,"r2_barplot.pdf"), r2.bars, base_width=8)



#############################################################################################
#################################### Generate scatterplot ###################################
#############################################################################################

full.spread %>% filter(repl == "rep0", sim == "sim100taxa") -> rep0.spread
            

plota <- rep0.spread %>% 
            ggplot(aes(x = homo_hyphy, y = ML_homo)) +
            geom_point() + 
            geom_abline(slope=1, intercept=0, color="blue") + 
            coord_cartesian(xlim=c(0,10), ylim=c(0,10))+ 
            xlab("LEISR") + ylab("R4S")

plotb <- rep0.spread %>% 
            ggplot(aes(x = gamma_hyphy, y = ML_gamma)) +
            geom_point() + 
            geom_abline(slope=1, intercept=0, color="blue") + 
            coord_cartesian(xlim=c(0,10), ylim=c(0,10)) +
            xlab("LEISR+G") + ylab("R4S+G")

            
plotc <- rep0.spread %>% 
            ggplot(aes(x = homo_hyphy, y = gamma_hyphy)) +
            geom_point() + 
            geom_abline(slope=1, intercept=0, color="blue") + 
            coord_cartesian(xlim=c(0,10), ylim=c(0,10)) +
            xlab("LEISR") + ylab("LEISR+G")
            
            
scatters <- plot_grid(plota, plotb, plotc, nrow=1, labels="auto")
save_plot(paste0(figpath,"replicate_100taxa_scatterplot.pdf"), scatters, base_width=8, base_height=2)


library(tidyverse)
theme_set(theme_classic() + theme(axis.title = element_text(size=11), axis.text = element_text(size=10)))

figpath = "../figures/"

r <- read.csv("r4s.csv")
l <- read.csv("leisr.csv")
tibble(leisr = l$MLE/mean(l$MLE), r4s = r$rate/mean(r$rate)) %>% 
    summarize(r = cor(leisr, r4s)^2) ## 0.9313825

tibble(leisr = l$MLE/mean(l$MLE), r4s = r$rate/mean(r$rate)) %>% 
    ggplot(aes(x = leisr, y=r4s)) +geom_point() +geom_abline(col="blue") + 
    xlab("LEISR") + ylab("R4S") -> p

ggsave(paste0(figpath,"hrh1_compare.pdf"), p)


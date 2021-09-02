library(tidyverse)

setwd("~/Github/carlson-betacov/")

bart.c <- (read.csv("CarlsonBartCitations.csv")[,-1] %>% rename(pred.bc = pred))
dart.c <- (read.csv("CarlsonDartCitations.csv")[,-1] %>% rename(pred.dc = pred))
bart.u <- (read.csv("CarlsonBartUncorrected.csv")[,-1] %>% rename(pred.bu = pred))
dart.u <- (read.csv("CarlsonDartUncorrected.csv")[,-1] %>% rename(pred.du = pred))

setwd("~/Github/Fresnel/Github/BatCSVs")

alb <- (read.csv("AlberyBats.csv")[,-1] %>% rename(pred.a = Count, host_species = Sp) %>% 
          mutate(host_species = gsub("_"," ",host_species)))
head(alb)

left_join(alb, bart.c) %>% left_join(dart.c) %>% left_join(bart.u) %>% left_join(dart.u) %>%
  select(host_species, betacov, pred.a, pred.bc, pred.dc, pred.bu, pred.du) -> compare

cor(compare[,2:7])          

library(WVPlots)

compare$betacov <- factor(compare$betacov)

colnames(compare)[4:7] <- c('BART-c', 'DART-c', 'BART-u', 'DART-u')
PairPlot(compare, 
         colnames(compare)[4:7], 
         "BART model comparison among four formulations", 
         group_var = "betacov",
         alpha = 0.35) + 
  theme_bw() +
  ggplot2::scale_color_manual(values = c('grey','red'))

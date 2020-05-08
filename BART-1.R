
setwd('~/Github/carlson-batcov')
set.seed(05082020)

library(BART)
library(fastDummies)
library(tidyverse)

varimp.pbart <- function(model, plot=TRUE) {
  
  names <- colnames(model$varcount)
  varimps <- colMeans(model$varcount/rowSums(model$varcount))
  
  var.df <- data.frame(names, varimps)
  
  var.df$names <- factor(var.df$names)
  var.df <- transform(var.df, names = reorder(names,
                                              -varimps))
  
  if(plot==TRUE){
    #g1 <- ggplot2::ggplot(var.df, aes(y=varimps, x=names)) +
    #  geom_bar(stat="identity", color="black") +
    #  theme(axis.text.x = element_text(angle = 45)) + 
    #  ylab("Relative importance") + theme_bluewhite()
    #print(g1)
    
    rel <- model$varcount/rowSums(model$varcount)
    colnames(rel) <- names
    
    rel %>% data.frame() %>% gather() %>%
      group_by(key) %>%
      summarise(mean = mean(value),
                sd = sd(value, na.rm = TRUE)) %>% 
      transform(Var = reorder(key, mean)) %>%
      ggplot(aes(x = Var, y = mean)) +
      geom_pointrange(aes(y = mean, x = Var, ymin = mean-sd, ymax = mean+sd),
                      color="#00AFDD") + 
      xlab(NULL) + ylab("Variable importance") + coord_flip() + 
      theme_bw() + theme(legend.position = "none",
                         axis.title.x = element_text(size=rel(1.3), vjust = -0.8),
                         axis.text.y = element_text(size=rel(1.4)),
                         plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                         panel.grid.minor = element_blank(),
                         panel.grid.major.x = element_line(color='grey',
                                                           linetype='dashed')) -> p
    
    print(p)
    
  }
  
  return(var.df)
  
}

read_csv('~/GitHub/cleanbats_betacov/clean data/BatCoV-assoc_compatible.csv') %>% filter(origin == 'Anthony') -> batcov
read_csv('~/GitHub/cleanbats_betacov/clean data/Han-BatTraits_compatible.csv') -> traits

# Add outcome variables 

batcov %>% mutate(betacov = as.numeric(virus_genus == 'Betacoronavirus'),
                  sarbecov = as.numeric(virus_subgenus == 'Sarbecovirus')) -> batcov

# Replace NA's for subgenus with 0 only if they're not betacoronaviruses
# (Addresses large number of unclassified betacoronaviruses)

batcov$sarbecov[is.na(batcov$sarbecov) & !(batcov$virus_genus == 'Betacoronavirus')] <- 0
batcov$sarbecov[is.na(batcov$sarbecov) & (batcov$virus_genus == 'Betacoronavirus')] <- -1

batcov %>% select(host_species, betacov, sarbecov) %>% unique -> batcov

# An internal function to summarize known associations

shorthand <- function(x) {
  if(sum(x==1)>0) {1} else {
  if(sum(x==-1)>0) {-1} else {
    0
  }
  }
}

batcov %>% group_by(host_species) %>% 
  summarize(betacov = max(betacov),
            sarbecov = shorthand(sarbecov)) -> batcov

# Create binomial names in the trait data

traits %>% mutate(host_species = paste(MSW05_Genus, MSW05_Species)) %>% 
  mutate() -> traits

# Add traits and associations

right_join(batcov, traits) %>% 
  mutate(betacov = replace_na(betacov, 0),
         sarbecov = replace_na(sarbecov, 0)) %>%
  mutate(sarbecov = na_if(sarbecov, -1)) -> batdf 

batdf <- data.frame(batdf)

# Turn categorical variable into columns

batdf %>% dummy_cols('ForStrat.Value') %>% 
  select(-ForStrat.Value) -> batdf

# Full model fails because there are too many NA's somewhere
# Drop any variable with > 50% NA's

varnums <- c(8:74)

varkeep <- 7 + unname(which(c(colSums(is.na(batdf[,8:74]))/nrow(batdf)) < 0.5))

# Updated table


model1 <- pbart(x.train = batdf[,varkeep],
                y.train = batdf[,'betacov'],
                sparse = FALSE,
                ntree = 200L,
                ndpost = 10000L)

model2 <- pbart(x.train = batdf[,varkeep],
                y.train = batdf[,'betacov'],
                sparse = TRUE,
                ntree = 200L,
                ndpost = 10000L)

varimp.pbart(model1)
varimp.pbart(model2)

batdf$pred1 <- colMeans(pnorm(model1$yhat.train))
batdf$pred2 <- colMeans(pnorm(model2$yhat.train))

# A density plot in the style of Becker

batdf %>% 
  ggplot(aes(pred2, 
             fill = factor(betacov), 
             colour = factor(betacov))) + 
  geom_density(alpha = 0.1)
 
# Top rankings 

batdf %>% 
  as_tibble() %>%
  filter(!(betacov == 1)) %>%
  select(host_species, pred2) %>%
  arrange(-pred2) %>% View()

# Get a 90% omission threshold

batdf %>% 
  as_tibble() %>%
  select(host_species, betacov, pred2) %>%
  data.frame() -> training

library(PresenceAbsence)

thresh <- optimal.thresholds(data.frame(training),
                             threshold = 10001,
                             opt.methods = 10,
                             req.sens = 0.9,
                             na.rm = TRUE)[1,2]

# How many new bats are above the threshold?

batdf %>% 
  as_tibble() %>%
  filter(!(betacov == 1)) %>%
  select(host_species, pred2) %>%
  arrange(-pred2) %>% 
  filter(pred2 > thresh) -> not.df
nrow(not.df)

# How's the AUC look

auc.roc.plot(data.frame(training))

# File for Greg

batdf %>% select(host_species,
                 betacov,
                 sarbecov,
                 pred2) %>% mutate(pred.bin = (pred2 > thresh),
                                   rank = rank(pred2)) %>%
  mutate(rank = (max(rank) - rank + 1)) %>%
  rename(pred = pred2) %>% as_tibble() -> bat.report

write.csv(bat.report, 'batcov-bart.csv')

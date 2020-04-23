
setwd('~/Github/carlson-batcov')
library(BART)
library(fastDummies)
library(tidyverse)

set.seed(69)

read_csv('BatCoV-assoc.csv') %>% filter(origin == 'Anthony') -> batcov
read_csv('Han-BatTraits.csv') -> traits

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

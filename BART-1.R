
setwd('C:/Users/cjcar/Desktop/BatCoV')
library(embarcadero)

batcov <- read_csv('BatCoV-assoc.csv')
traits <- read_csv('Han-BatTraits.csv')

# Add outcome variables 

batcov %>% mutate(betacov = as.numeric(virus_genus == 'Betacoronavirus'),
                  sarbecov = as.numeric(virus_subgenus == 'Sarbecovirus')) -> batcov

# Replace NA's for subgenus with 0 only if they're not betacoronaviruses
# (Addresses large number of unclassified betacoronaviruses)

batcov$sarbecov[is.na(batcov$sarbecov) & !(batcov$virus_genus == 'Betacoronavirus')] <- 0
batcov$sarbecov[is.na(batcov$sarbecov) & (batcov$virus_genus == 'Betacoronavirus')] <- -1

batcov %>% select(host_species, betacov, sarbecov) %>% unique -> batcov

# Create binomial names in the trait data

traits %>% mutate(host_species = paste(MSW05_Genus, MSW05_Species)) %>% 
  mutate() -> traits

# Add traits and associations

right_join(batcov, traits) %>% 
  mutate(betacov = replace_na(betacov, 0),
         sarbecov = replace_na(sarbecov, 0)) %>%
  mutate(betacov = na_if(betacov, -1),
         sarbecov = na_if(sarbecov, -1)) -> batdf 

batdf <- data.frame(batdf)

# Turn categorical variable into columns

library(fastDummies)

batdf %>% dummy_cols('ForStrat.Value') %>% 
  select(-ForStrat.Value) -> batdf

# Full model fails because there are too many NA's somewhere
# Drop any variable with > 50% NA's

varnums <- c(8:74)

varkeep <- 7 + unname(which(c(colSums(is.na(batdf[,8:74]))/nrow(batdf)) < 0.9))

# Updated table

library(BART)

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

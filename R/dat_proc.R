library(tidyverse)
library(readr)
# library(ASCI)
devtools::load_all('../ASCI')


######
# pre-processing of Susie's raw data

alg <- read.csv('data/algae.bug.data.092317.csv', stringsAsFactors = F)
env <- read.csv('data/algae.site.data.092317.csv', stringsAsFactors = F)

algprc <- getids(alg)
envprc <- getids(env)

# check site matches
chkinp(algprc, envprc)
rmv <- chkinp(algprc, envprc, getval = T)
algprc <- algprc %>%
  filter(!SampleID %in% rmv)
envprc <- envprc %>%
  filter(!SampleID %in% rmv)

# check taxonomy
chkinp(algprc, envprc)
rmv <- chkinp(algprc, envprc, getval = T)
algprc <- algprc %>%
  filter(!FinalID %in% rmv)

# check if both sba, diatoms are present
chkinp(algprc, envprc)
rmv <- chkinp(algprc, envprc, getval = T)
algprc <- algprc %>%
  filter(!SampleID %in% rmv)
envprc <- envprc %>%
  filter(!SampleID %in% rmv)

# check diatom abundance data
chkinp(algprc, envprc)
rmv <- chkinp(algprc, envprc, getval = T)
algprc <- algprc %>%
  filter(!SampleID %in% rmv)
envprc <- envprc %>%
  filter(!SampleID %in% rmv)

# check NA values in environmental predictors
chkinp(algprc, envprc)
rmv <- chkinp(algprc, envprc, getval = T) %>% 
  .$SampleID %>% 
  unique
algprc <- algprc %>%
  filter(!SampleID %in% rmv)
envprc <- envprc %>%
  filter(!SampleID %in% rmv)

# final
dat <- chkinp(algprc, envprc)

# get asci object and save
results <- ASCI(dat$taxa, dat$site)
save(results, file = 'data/results.RData', compress = 'xz')

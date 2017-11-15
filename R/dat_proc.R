library(tidyverse)
library(readr)
# library(ASCI)
devtools::load_all('../ASCI')

alg <- read.csv('data/algae.bug.data.092317.csv', stringsAsFactors = F)
env <- read.csv('data/algae.site.data.092317.csv', stringsAsFactors = F)



algprc <- getids(alg)

rmv <- chkinp(algprc, env, getval = T)
algprc <- algprc %>%
  filter(!SampleID %in% rmv)
rmv <- chkinp(algprc, env, getval = T)
algprc <- algprc %>%
  filter(!FinalID %in% rmv)
rmv <- chkinp(algprc, env, getval = T)
algprc <- algprc %>% 
  filter(!SampleID %in% rmv)




ASCI(algprc, env)

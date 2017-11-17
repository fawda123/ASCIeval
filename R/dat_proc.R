library(tidyverse)
library(readr)
library(tidyverse)
library(forcats)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(readxl)
library(lubridate)

library(ASCI)
# devtools::load_all('../ASCI')

######
# pre-processing of Susie's raw data

# alg <- read.csv('ignore/algae.bug.data.092317.csv', stringsAsFactors = F)
# env <- read.csv('ignore/algae.site.data.092317.csv', stringsAsFactors = F)
# 
# algprc <- getids(alg)
# envprc <- getids(env)
# 
# # check site matches
# chkinp(algprc, envprc)
# rmv <- chkinp(algprc, envprc, getval = T)
# algprc <- algprc %>%
#   filter(!SampleID %in% rmv)
# envprc <- envprc %>%
#   filter(!SampleID %in% rmv)
# 
# # check taxonomy
# chkinp(algprc, envprc)
# rmv <- chkinp(algprc, envprc, getval = T)
# algprc <- algprc %>%
#   filter(!FinalID %in% rmv)
# 
# # check if both sba, diatoms are present
# chkinp(algprc, envprc)
# rmv <- chkinp(algprc, envprc, getval = T)
# algprc <- algprc %>%
#   filter(!SampleID %in% rmv)
# envprc <- envprc %>%
#   filter(!SampleID %in% rmv)
# 
# # check diatom abundance data
# chkinp(algprc, envprc)
# rmv <- chkinp(algprc, envprc, getval = T)
# algprc <- algprc %>%
#   filter(!SampleID %in% rmv)
# envprc <- envprc %>%
#   filter(!SampleID %in% rmv)
# 
# # check NA values in environmental predictors
# chkinp(algprc, envprc)
# rmv <- chkinp(algprc, envprc, getval = T) %>% 
#   .$SampleID %>% 
#   unique
# algprc <- algprc %>%
#   filter(!SampleID %in% rmv)
# envprc <- envprc %>%
#   filter(!SampleID %in% rmv)
# 
# # final and save
# dat <- chkinp(algprc, envprc)
# taxain <- dat$taxa
# sitein <- dat$site
# 
# save(taxain, file = 'data/taxain.RData')
# save(sitein, file = 'data/sitein.RData')

# get asci object and save
data(taxain)
data(sitein)
pkgdat <- ASCI(taxain, sitein)
save(pkgdat, file = 'data/pkgdat.RData', compress = 'xz')

######
# site meta - ref/int/str, cal/val

##
# import environmental data and scores for OE, MMI, clip all by SMC boundaries

# all data
fls <- list.files('ignore', pattern = 'Scores.*\\.csv$', full.names = TRUE) %>% 
  tibble(fl = .) %>% 
  mutate(
    data = map(fl, read.csv, stringsAsFactors = FALSE)
  )

# OE data, filter by sites in latlon
sitcat <- fls %>% 
  unnest %>%
  dplyr::select(X, Type) %>% 
  filter(!Type %in% 'notrecent') %>%
  unique %>% 
  rename(SampleID = X) %>% 
  mutate(SampleID = gsub('(_[0-9]+)\\.([0-9]+)\\.([0-9]+_)', '\\1/\\2/\\3', SampleID))  

save(sitcat, file = 'data/sitcat.RData', compress = 'xz')

######
# original data for checking

# all data
fls <- list.files('legacy', pattern = '\\.csv$', full.names = TRUE) %>% 
  tibble(fl = .) %>% 
  mutate(
    data = map(fl, read.csv, stringsAsFactors = FALSE)
  )

# OE data, filter by sites in latlon
oedat <- fls %>% 
  filter(grepl('Scores', fl)) %>% 
  unnest %>% 
  rename(
    SampleID = X,
    scr = OE.scores.OoverE
  ) %>% 
  mutate(
    grp = gsub('legacy/Scores\\.OE\\.|\\.csv*$', '', fl),
    grp = gsub('diatom$', 'diatoms', grp),
    ind = 'oe'
  ) %>% 
  dplyr::select(SampleID, ind, grp, scr)

# MMI data
mmdat <- fls %>% 
  filter(grepl('results', fl)) %>% 
  unnest %>% 
  dplyr::select(X, matches('pMMI')) %>% 
  gather('grp', 'scr', -X) %>% 
  na.omit %>% 
  mutate(
    grp = gsub('hybri', 'hybrid', grp), 
    grp = gsub('\\.pMMI', '', grp), 
    grp = gsub('diatom$', 'diatoms', grp),
    ind = 'mmi'
  ) %>% 
  rename(
    SampleID = X
  ) %>% 
  dplyr::select(SampleID, ind, grp, scr)

# combine mmdat, oedat
orgdat <- rbind(mmdat, oedat) %>% 
  mutate(dat = 'org')

save(orgdat, file = 'data/orgdat.RData', compress = 'xz')

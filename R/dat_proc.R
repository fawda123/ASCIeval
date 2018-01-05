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
library(sf)
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

######
# psa regions 

prstr <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
psa <- st_read('S:/Spatial_Data/RCMP_needs editting/Inputs/PSA6_090111/PSA6_2011.shp') %>% 
  st_transform(prstr)
save(psa, file = 'data/psa.RData', compress = 'xz')

######
# site meta - ref/int/str, cal/val, psa region

##
# import environmental data and scores for OE, MMI, clip all by SMC boundaries

# all data
fls <- list.files('ignore', pattern = 'Scores.*\\.csv$', full.names = TRUE) %>% 
  tibble(fl = .) %>% 
  mutate(
    data = map(fl, read.csv, stringsAsFactors = FALSE)
  )

# site category - refernece (cal, val), stressed, intermediate, and notrecent (repeats)
sitcat <- fls %>% 
  unnest %>%
  dplyr::select(X, Type) %>%
  unique %>% 
  rename(
    SampleID = X,
    cls = Type
    ) %>% 
  mutate(SampleID = gsub('(_[0-9]+)\\.([0-9]+)\\.([0-9]+_)', '\\1/\\2/\\3', SampleID))  

# intersect site locations with psa regions
sitmet <- sitein %>% 
  dplyr::select(SampleID, New_Long, New_Lat) %>% 
  st_as_sf(coords = c('New_Long', 'New_Lat')) %>% 
  st_set_crs(prstr) %>% 
  st_intersection(psa)

# remove geometry
sitmet_nogeo <- sitmet %>% 
  dplyr::select(SampleID, PSA_REGION) %>% 
  st_set_geometry(NULL) %>% 
  rename(psa = PSA_REGION)

sitcat <- sitcat %>% 
  left_join(sitmet_nogeo, by = 'SampleID')

save(sitcat, file = 'data/sitcat.RData', compress = 'xz')

##
# other algae ibi scores

# old algal ibi divided by median of ref calibration scores to standardize
thrsh <- tibble(
  ind = c('H20','D18','S2'),
  thrsh = c(75, 79, 69)
)

# process
aldat <- read.csv('ignore/tblAlgaeIBI.csv', stringsAsFactors = FALSE) %>% 
  dplyr::select(StationCode, SampleDate, Replicate, S2, D18, H20) %>% 
  mutate(
    SampleDate = mdy_hms(SampleDate, tz = 'Pacific/Pitcairn'),
    SampleMonth = month(SampleDate),
    SampleDay = day(SampleDate), 
    SampleYear = year(SampleDate)
  ) %>% 
  dplyr::select(-SampleDate) %>% 
  unite('SampleDate', SampleMonth, SampleDay, SampleYear, sep ='/') %>% 
  unite('SampleID', StationCode, SampleDate, Replicate, sep = '_') %>% 
  filter(!duplicated(SampleID)) %>% # one duplicated....
  gather('ind', 'scr', S2:H20) %>%
  left_join(thrsh, by = 'ind') %>% 
  mutate(
    scr = scr / thrsh
    ) %>% 
  dplyr::select(-thrsh)

# save
save(aldat, file = 'data/aldat.RData', compress = 'xz')

######
# save a copy of taxain processed to genus or higher

data(taxain)

# this creates a lookup table from pmmilkup$traits in asci package 
x <- pmmilkup$traits %>%
  dplyr::select(Phylum, Class, Order, Family, Genus1, FinalIDassigned) %>% 
  mutate_if(is.factor, as.character) %>% 
  rowwise() %>% 
  mutate(Genus2 = do({
      
      # use genus1 if not empty
      if(Genus1 != '')
        return(Genus1)
      
      chk <- c(Phylum, Class, Order, Family, Genus1)
      
      # get first of FinalIDassigned if all higher taxa cats are empty
      if(all(chk == '')){
        
        out <- strsplit(FinalIDassigned, ' ')[[1]][1]
        
      # get last non-empty cat  
      } else {
        
        out <- chk[chk != ''] %>% rev %>% .[1]
        
      }
      
      return(out)
      
      })
    ) %>% 
  dplyr::select(FinalIDassigned, Genus2) %>% 
  rename(
    FinalID = FinalIDassigned, 
    Genus = Genus2
    )

# join taxain with key, then use strsplit
taxain_gen <- taxain %>% 
  left_join(x, by = 'FinalID') %>% 
  mutate(Genus = ifelse(is.na(Genus), gsub(' .*$', '', FinalID), Genus)) %>% 
  dplyr::select(-FinalID) %>% 
  rename(FinalID = Genus)

# save
save(taxain_gen, file = 'data/taxain_gen.RData', compress = 'xz')

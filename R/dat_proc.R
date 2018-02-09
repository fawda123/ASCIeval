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
# taxain_spp <- dat$taxa
# sitein <- dat$site
# 
# save(taxain_spp, file = 'data/taxain_spp.RData')
# save(sitein, file = 'data/sitein.RData')

# get asci object and save
data(taxain_spp)
data(sitein)
pkgdat <- ASCI(taxain_spp, sitein)
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
# save a copy of taxain_spp processed to genus or higher

data(taxain_spp)
data(pmmilkup_spp)

# this creates a lookup table from pmmilkup$traits in asci package 
x <- pmmilkup_spp$traits %>%
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

# join taxain_spp with key, then use strsplit
taxain_gen <- taxain_spp %>% 
  left_join(x, by = 'FinalID') %>% 
  mutate(Genus = ifelse(is.na(Genus), gsub(' .*$', '', FinalID), Genus)) %>% 
  dplyr::select(-FinalID) %>% 
  rename(FinalID = Genus) %>% 
  group_by(SampleID, StationCode, SampleDate, Replicate, SampleTypeCode, FinalID) %>% 
  summarise(
    Result = sum(Result, na.rm = T), 
    Result = ifelse(Result == 0, NA, Result),
    BAResult = sum(BAResult, na.rm = T),
    BAResult = ifelse(BAResult == 0, NA, BAResult)
  ) %>% 
  dplyr::select(SampleID, StationCode, SampleDate, Replicate, SampleTypeCode, BAResult, Result, FinalID)

# save
save(taxain_gen, file = 'data/taxain_gen.RData', compress = 'xz')

######
# create genus level demo algae taxonomy data

data(demo_algae_tax_spp)
data(pmmilkup_spp)

# this creates a lookup table from pmmilkup$traits in asci package 
x <- pmmilkup_spp$traits %>%
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

# join taxain_spp with key, then use strsplit
demo_algae_tax_gen <- demo_algae_tax_spp %>% 
  left_join(x, by = 'FinalID') %>% 
  mutate(Genus = ifelse(is.na(Genus), gsub(' .*$', '', FinalID), Genus)) %>% 
  dplyr::select(-FinalID) %>% 
  rename(FinalID = Genus) %>% 
  group_by(StationCode, SampleDate, Replicate, CollectionMethodCode, SampleTypeCode, FinalID) %>% 
  summarise(
    Result = sum(Result, na.rm = T), 
    Result = ifelse(Result == 0, NA, Result),
    BAResult = sum(BAResult, na.rm = T),
    BAResult = ifelse(BAResult == 0, NA, BAResult)
  ) %>% 
  dplyr::select(StationCode, SampleDate, Replicate, CollectionMethodCode, SampleTypeCode, BAResult, Result, FinalID)

# save
save(demo_algae_tax_gen, file = 'data/demo_algae_tax_gen.RData', compress = 'xz')

######
# create a genus level trait table, must export to ASCI package genus branch

# species traits table with complete genus column
data(pmmilkup_spp)
data(taxain_spp)
traits_gen <- trt_gen_fun(taxain_spp, pmmilkup_spp$traits, thrsh = 0.75)

# save(traits_gen, file = 'C:/Users/Marcus.SCCWRP2K/Desktop/traits_gen.RData', compress = 'xz')

######
# indicator species analysis at genus level
library(tidyverse)
library(indicspecies)
library(ASCI)

data(taxain_gen)
data(sitcat)

##
# groups for indic
grpin <- sitcat %>% 
  filter(cls %in% c('int', 'rc', 'str')) %>% 
  mutate(
    cls = fct_recode(cls, Intermediate = 'int',  Reference = 'rc', Stressed = 'str'),
    cls = as.character(cls)
  ) %>% 
  dplyr::select(-psa) %>% 
  spread(SampleID, cls) %>% 
  unlist

##
# filter taxonomy data by sites in dev dataset

# format abundance ests
abuin <- taxain_gen %>% 
  dplyr::select(SampleID, FinalID, BAResult, Result) %>% 
  mutate(
    abund = as.numeric(pmax(BAResult, Result, na.rm = T))
  ) %>% # combine abundance 
  dplyr::select(-BAResult, -Result) %>% 
  group_by(SampleID, FinalID) %>% 
  summarise(abund = sum(abund, na.rm = T)) %>%
  ungroup

# make wide format, location to rowid
abuin <- abuin %>% 
  spread(FinalID, abund, fill = 0) %>% 
  data.frame(stringsAsFactors = F) %>% 
  remove_rownames %>% 
  column_to_rownames('SampleID')

##
# make sure sites match between grpin and abuin, order is the same
grpin <- grpin[names(grpin) %in% rownames(abuin)]
abuin <- abuin %>% 
  .[rownames(.) %in% names(grpin), ] %>% 
  .[match(rownames(.), names(grpin)), ]

##
# run indicator analysis
res <- multipatt(abuin, grpin, duleg = T, control = how(nperm = 999))

# the A value describes the affinity of a species for each group
# this sums to one for each species across the groups
# it indicates which group a species is most likely to belong
# the B value describes the proportion of sites within a group that have a species, if found in the group
# the final stat value for indic analysis considers A and B
# multipatt returns a pvalue for each species
# the species indicators can be further selected by filtering with an A or B threshold
# the coverage function describes the proportion of sites in each group where one or
# indicators are found for a given threshold, with A = 0 being the default
# using a more conservative, larger A will remove indicator species but this can be done
# up to a threshold until coverage begins to decline

# table of indicator species by pvalue
indic <- res$sign %>%
  rownames_to_column('FinalID') %>% 
  filter(p.value < 0.05) %>% 
  mutate(designation = factor(index, 
                              levels = c('1', '2', '3'), 
                              labels = c('Intermediate', 'Reference', 'Stressed'))
         
  ) %>% 
  select(FinalID, designation, stat, p.value)

# all a values
aval <- res$A %>% 
  data.frame %>% 
  rownames_to_column('FinalID') %>% 
  gather('designation', 'aval', -FinalID)

# find maximum a values without decreasing coverage
chks <- seq(0, 1, by = 0.005)
toplo <- sapply(chks, function(At) coverage(x = abuin, y = res, At = At)) %>% 
  t %>% 
  data.frame %>% 
  mutate(chks = chks) %>% 
  gather('designation', 'cvg', -chks)
athr <- toplo %>% 
  group_by(designation) %>%
  filter(cvg == max(cvg)) %>% 
  filter(chks == max(chks)) %>% 
  select(-cvg) %>% 
  ungroup
# ggplot(toplo, aes(x = chks, y = cvg, group = designation, colour = designation)) +
#   geom_vline(data = athr, aes(xintercept = chks, colour = designation), linetype = 'dashed') +
#   scale_x_continuous('At') + 
#   geom_line() +
#   theme_bw()

# join indic with a val, a threshold
indic <- indic %>% 
  left_join(aval, by = c('FinalID', 'designation')) %>% 
  left_join(athr, by = c('designation')) %>% 
  filter(aval >= chks)

# final
indicators <- indic %>%
  dplyr::select(FinalID, designation) %>%
  rename(FinalIDassigned = FinalID)

# save(indicators, file = 'C:/Users/Marcus.SCCWRP2K/Desktop/indicators.RData', compress = 'xz')

#####
# estimate raw metrics at genus and species level for algal data

library(tidyverse)
library(reshape2)
library(ASCI)
library(vegan)
source('R/funcs.R')

data(taxain_gen)
data(taxain_spp)
data(pmmilkup_gen)
data(pmmilkup_spp)
data(sitein)
data(STE)

indicators_gen <- pmmilkup_gen$indicators
traits_gen <- pmmilkup_gen$traits
indicators_spp <- pmmilkup_spp$indicators
traits_spp <- pmmilkup_spp$traits

rawmet_gen <- pmmifun_hrd(taxain_gen, sitein, indicators_gen, traits_gen)
rawmet_spp <- pmmifun_hrd(taxain_spp, sitein, indicators_spp, traits_spp)

save(rawmet_gen, file = 'data/rawmet_gen.RData', compress = 'xz')
save(rawmet_spp, file = 'data/rawmet_spp.RData', compress = 'xz')

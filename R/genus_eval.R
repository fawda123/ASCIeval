library(tidyverse)
library(reshape2)
library(ASCI)
library(vegan)
library(foreach)
library(doParallel)
source('R/funcs.R')

data(taxain_gen)
data(sitein)
data(STE)
data(pmmilkup_spp)
data(taxain_spp)
data(rawmet_spp)

# metrics to evaluate
toeval <- pmmilkup_spp$quants %>% .$Metric %>% unique %>% sort

# proportion thresholds to evaluate
prps <- seq(0.25, 1, by = 0.02)

# setup parallel backend
cl<-makeCluster(6)
registerDoParallel(cl)
strt<-Sys.time()

# get metric estimates for each proportion threshold
out_ls <- foreach(i = prps, .packages = c('ASCI', 'tidyverse', 'reshape2', 'vegan')) %dopar% {

  source('R/funcs.R')
  sink('log.txt')
  cat(i, which(i == prps), 'of', length(prps), '\n')
  print(Sys.time()-strt)
  sink()
  
  load(file = 'data/pmmilkup_spp.RData')
  load(file = 'data/pmmilkup_gen.RData')
  indicators_gen <- pmmilkup_gen$indicators
  
  traits_gen <- trt_gen_fun(taxain_spp, pmmilkup_spp$traits, thrsh = i)
  
  rawmet_gen <- pmmifun_hrd(taxain_gen, sitein, indicators_gen, traits_gen, mets = toeval)
  
  return(rawmet_gen)
  
}

# lookup table for getting metrics specific to each taxa
tojn <- pmmilkup_spp$quants %>% 
  select(Metric, Assemblage) %>% 
  rename(
    name1 = Assemblage,
    met = Metric
    )

# species trait metrics for comparison
met_spp <- rawmet_spp %>% 
  select(-res) %>% 
  filter(met %in% toeval)

# output formatted for select metrics by taxa
# joined w/ 'actual' species metrics
ests <- out_ls %>%
  enframe %>%
  mutate(
    prps = prps
    ) %>%
  unnest %>%
  filter(!met %in% c('richness', 'shannon', 'simpson')) %>% 
  left_join(tojn, ., by = c('met', 'name1')) %>% 
  select(-name) %>%
  rename(name = name1) %>% 
  left_join(met_spp, by = c('met', 'name', 'SampleID')) %>% 
  rename(
    val_gen = val.x, 
    val_spp = val.y
  )
# percentages to sum for genus, species metric comparison
prchk <- seq(0, 0.5, by = 0.05)

# summarize difference between gen and spp mets
sms <- ests %>% 
  group_by(met, name, prps) %>% 
  nest %>% 
  mutate(
    dff = map(data, function(x){
      browser()
      
      lm()
      prdf <- x %>% 
        mutate(
          prdf = abs(val_gen - val_spp)/spp_rng
        ) %>% 
        .$prdf 
      
      out <- sapply(prchk, function(x) sum(prdf <= x))
      names(out) <- prchk
      out <- data.frame(cnt = out) %>% 
        rownames_to_column('wthn')
        
      return(out)
      
      
    })
  ) %>% 
  select(-data) %>% 
  unnest

ggplot(sms, aes(x = prps, y = cnt/ 2300, group = wthn, colour = wthn)) +
  facet_wrap(name ~ met) +
  theme_minimal() + 
  theme(
    panel.background = element_rect(fill = 'lightgrey'),
    panel.border = element_rect(fill = NA), 
    panel.grid.minor = element_blank()
    ) + 
  geom_line() + 
  scale_colour_brewer(palette = 'Spectral') + 
  scale_x_continuous('Minimum threshold individuals in a trait') +
  scale_y_continuous('Number of sites') + 
  geom_hline(yintercept = 0.75, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed')
  
  
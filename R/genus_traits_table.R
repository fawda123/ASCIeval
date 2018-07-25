data(taxain_spp)
traits_in <- read.csv('ignore/algae.traits.data.07242018.csv', stringsAsFactors = F)

# add a complete genus column to the traits lookup table
spp_gen_trt <- traits_in %>%
  dplyr::select(Phylum, Class, Order, Family, Genus, FinalIDassigned, everything()) %>% 
  dplyr::select(-AlgaeList, -FinalID, -BCGfinal) %>% 
  mutate_if(is.factor, as.character) %>% 
  rowwise() %>% 
  mutate(Genus = do({

    # use Genus if not empty
    if(!is.na(Genus))
      return(Genus)
    
    chk <- c(Phylum, Class, Order, Family, Genus)
    
    # get first of FinalIDassigned if all higher taxa cats are empty
    if(all(is.na(chk))){
      
      out <- strsplit(FinalIDassigned, ' ')[[1]][1]
      
      # get last non-empty cat  
    } else {
      
      out <- chk[!is.na(chk)] %>% rev %>% .[1]
      
    }
    
    return(out)
    
  })
  ) %>% 
  dplyr::select(-Phylum, -Class, -Order, -Family) %>% 
  unique %>% 
  gather('trait', 'est', -Genus, -FinalIDassigned) %>% 
  rename(FinalID = FinalIDassigned)

# join species level data with genus traits to summarize counts by data dist
# genus from unmatched finalid in taxain_spp added to genus in spp_gen_trt
to_summ <- taxain_spp %>% 
  left_join(spp_gen_trt, by = 'FinalID') %>% 
  mutate(
    Genus = ifelse(is.na(Genus), gsub(' .*$', '', FinalID), Genus),
    est = ifelse(est %in% '', NA, est)
  )

thrsh <- c(0.5, 0.75, 1)

for(i in thrsh){
  
  # summarize trait counts by data dist, get those with greater than or equal to thrsh, spread
  traits_gen <- to_summ %>% 
    group_by(Genus, trait, est) %>% 
    summarise(n = length(est)) %>% 
    group_by(Genus, trait) %>% 
    mutate(prp = n / sum(n)) %>% 
    group_by(Genus, trait) %>% 
    filter(prp == max(prp)) %>% # sometimes there are multiple trait categories above thrsh, pick max (most likely)
    filter(!duplicated(prp)) %>% # sometimes the maximum proportions are equal, just pick the first one
    filter(prp >= i) %>% # filter proportions that are above the hreshold
    ungroup %>% 
    dplyr::select(-prp, -n) %>% 
    spread(trait, est) %>% 
    rename(FinalIDassigned = Genus)

  flout <- i %>% 
    gsub('\\.', '_', .) %>% 
    paste0('ignore/traits', ., '.csv')
  
  write.csv(traits_gen, flout, row.names = F, na = '')
  
}
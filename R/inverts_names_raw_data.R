
# Check on field identification for invertebrates to see which names were used and how often they occurred



## load libaries

library( tidyverse )


## read files

# Data
d <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data.csv" )


barnacles <- c( "Semibalanus cariosus", "Pollicipes polymerus", "Chthamalus dalli",  "Barnacle", "Balanus sp.", "Balanus glandula")
anemones <- c( "Anemone", "Anthopleura elegantissima", "Anthopleura sp.", "Anthopleura xanthogrammica" )
mussels <- c( "Mytilus sp.", "Mytilus trossulus" )

invert_list <- list( barnacles, anemones, mussels )

result <- list() 
for(i in 1:length(invert_list)){
  result[[i]] <- d %>% separate( UID, c("Site", "Beach", "Zone", "Year", "Quadrat") ) %>% 
    filter( Year %in% 2012:2019) %>% 
    filter( Site != "Meay Channel") %>% 
    filter( Taxon %in% invert_list[[i]] ) %>% 
    group_by(Taxon) %>% 
    summarize( occurrences = length(Abundance), mean_cover_when_present = mean(Abundance) ) %>% 
    mutate( proportional_occurrence = occurrences/sum(occurrences), 
            proportional_cover = mean_cover_when_present/sum(mean_cover_when_present) ) %>% 
    arrange( -proportional_occurrence )
}
write.csv( do.call(rbind, result), "R/output/invert_field_names_representation.csv")

invert_vector <- do.call(c, invert_list)
site_taxon <- d %>% 
  separate( UID, c("Site", "Beach", "Zone", "Year", "Quadrat") ) %>% 
  filter( Year %in% 2012:2019) %>% 
  filter( Site != "Meay") %>% 
  filter( Taxon %in%  invert_vector) %>% 
  group_by(Site, Taxon) %>% 
  summarize( occurrences = length(Abundance), mean_cover_when_present = mean(Abundance) )
write.csv( site_taxon, "R/output/invert_field_names_representation_sites.csv")

# barnacles by site and zone
barnacle_zone <- d %>% 
  separate( UID, c("Site", "Beach", "Zone", "Year", "Quadrat") ) %>% 
  mutate( Zone = factor(Zone, levels = c("HIGH","MID","LOW"))) %>% 
  filter( Year %in% 2012:2019) %>% 
  filter( Site != "Meay") %>% 
  filter( Taxon %in%  barnacles) %>% 
  group_by(Site, Zone, Taxon) %>% 
  summarize( occurrences = length(Abundance), mean_cover_when_present = mean(Abundance) )
write.csv( barnacle_zone, "R/output/invert_field_names_representation_barnacle_zone.csv")



# Look at breakdowns just in 2012 and 2014 when we were more careful to record species
result <- list() 
for(i in 1:length(invert_list)){
  result[[i]] <- d %>% separate( UID, c("Site", "Beach", "Zone", "Year", "Quadrat") ) %>% 
    filter( Year %in% c(2012,2014)) %>% 
    filter( Site != "Meay Channel") %>% 
    filter( Taxon %in% invert_list[[i]] ) %>% 
    group_by(Taxon) %>% 
    summarize( occurrences = length(Abundance), mean_cover_when_present = mean(Abundance) ) %>% 
    mutate( proportional_occurrence = occurrences/sum(occurrences), 
            proportional_cover = mean_cover_when_present/sum(mean_cover_when_present) ) %>% 
    arrange( -proportional_occurrence )
}
write.csv( do.call(rbind, result), "R/output/invert_field_names_representation_2012.csv")


# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space


# load libraries
library( tidyverse )
library( here )


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )



## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
# muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]

# # no NA values allowed, so we need to remove these from the dataset
# rem <- unique( which(is.na(muse$Shore_height_cm))  )
# # get rid of row 2 in all data structures
# muse   <- muse[-rem,]


muse <- droplevels(muse)


# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to sessile organisms (mobile data not very reliable)
d <- dm %>% 
  filter( motile_sessile != "motile" )
  # filter( non.alga.flag =="Algae" )

# set factor for zone
d$Zone <- factor( d$Zone, levels=c('LOW','MID','HIGH'))

# lump some taxa
dupdate <- d %>% filter( Taxon != taxon_revised )
d$taxon_revised <- gsub( "Mastocarpus.*", "Mastocarpus sp.", d$taxon_revised)
d$taxon_revised <- gsub( "Hildenbrandia.*", "Hildenbrandia sp.", d$taxon_revised)


# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  # mutate( taxon = gsub(" ",".",taxon_lumped3) ) %>% 
  mutate( taxon = gsub(" ",".",taxon_revised) ) %>% 
  group_by( UID, Year, Site, Zone, taxon, non.alga.flag ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) 

# spread and gather to fill with zeros
dzeros <- d.simple %>% 
  select( -non.alga.flag ) %>% 
  group_by(UID, Year, Site, Zone ) %>% 
  spread( taxon, Abundance, fill=0 ) %>% 
  gather( taxon, Abundance, -UID, -Year, -Site, -Zone  )

# look at a single year
drank <- dzeros %>% 
  filter( Year == 2019 ) %>% 
  group_by( Site, Zone, taxon ) %>% 
  summarize( Abundance=mean(Abundance,na.rm=T)) %>% 
  mutate( rank = rank( -Abundance, ties.method = "first" )) %>% 
  group_by( Site, Zone ) %>% 
  mutate( total=sum(Abundance), prop=Abundance/total ) %>% 
  arrange(Site, Zone, -prop ) %>% 
  mutate( csum = cumsum(prop), total = round(total,1), 
          prop = round(prop,2), csum = round(csum,2) )
dtop <- drank %>% 
  filter( rank <=15 )
dtop$taxon <- gsub( "[.]"," ", dtop$taxon )
dtop$taxon <- gsub( "sp ","sp.", dtop$taxon )
write_csv( dtop,"R Code and Analysis/output from r/ranks_2019.csv")
formattable::formattable( dtop, align = c("l","l","l")  )
# 
# drankt <- drank %>% 
#   group_by( Site, Zone ) %>% 
#   summarize( taxa = list(unique(taxon[order(rank)])) )
# # drank3t$taxa <- do.call( rbind, lapply( drank3t$taxa, function(z) paste(unlist(z),sep='',collapse=", ") ) )
# widetax <- data.frame( do.call(rbind,drankt$taxa) )
# colnames(widetax) <- paste0("rank",1:15)
# drankt <- bind_cols( select(drankt,-taxa), widetax )
# drankt$taxa <- gsub( "[.]"," ", drankt$taxa )
# # drank3t$taxa <- gsub( " sp "," sp.", drank3t$taxa )
# 
# write_csv( drankt,"R Code and Analysis/output from r/ranks_2019.csv")
# formattable::formattable( drankt, align = c("l","l","l")  )



# do all years to help make datasheets for future use
dking <- d.simple %>% ungroup() %>% 
  select( taxon, non.alga.flag ) %>% distinct()
dall <- left_join( dzeros, dking )
drankall <- dall %>% 
  # filter( Year == 2019 ) %>% 
  group_by( Site, Zone, non.alga.flag, taxon ) %>% 
  summarize( Abundance=mean(Abundance,na.rm=T)) %>% 
  mutate( rank = rank( -Abundance, ties.method = "first" )) %>% 
  group_by( Site, Zone ) %>% 
  mutate( total=sum(Abundance), prop=Abundance/total ) %>% 
  arrange(Site, Zone, non.alga.flag, -prop ) %>% 
  mutate( csum = cumsum(prop), total = round(total,1), 
          prop = round(prop,2), csum = round(csum,2) ) %>% 
  filter( Abundance > 0 )
write_csv( drankall,"R Code and Analysis/output from r/ranks_all.csv")

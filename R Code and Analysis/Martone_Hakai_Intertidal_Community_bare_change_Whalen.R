# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen (adapted in part from Sam Starko's diversity script)
# created 26 Feb 2019

# This script is SS's attempt at running community analyses 
# This script uses data on ALGAE ONLY from previous scripts (no animals, no substrate types etc)

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(vegan)
library(bayou)
library(viridis)
#run summarize functions below before loading Rmisc. 
#calling Rmisc seems to block summarize from doing its job.
# library(Rmisc)




## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read.csv( "Data/R Code/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]

# remove taxa that are not coutned towards subtratum cover (i.e. mobile invertebrates)
# make it easier by replacing NA values for substratum
ds <- ad
ds$motile_sessile[ is.na(ds$motile_sessile) ] <- "Substratum"
ds <- ds[ ds$motile_sessile!="motile", ]


# remove bare space? Not yet

d <- ds



# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which wil reduce the size of the dataset a bit
d.simple <- d %>%
  group_by( UID, taxon_lumped, non.alga.flag, motile_sessile, Kelp.Fucoid.Turf, Littler.defined ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))



# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
muse  <- am
drestrict <- d.simple[ d.simple$UID %in% muse$UID , ] 
dm <- left_join( drestrict, am )



# look at bar rock cover, barnacles, mussels across sites and zones
sub <- dm %>%
  group_by( UID, Date, Quadrat, Meter.point, Site, Zone, Year, Shore_height_cm, taxon_lumped ) %>%
  filter( taxon_lumped %in% c('Bare rock','Barnacles','Mytilus sp.')) %>%
  summarize( Abundance=sum(Abundance,na.rm=T) )
# get back abundances for other UIDs in west and fill with zero
all.sub  <- sub %>% ungroup() %>% tidyr::expand( UID, taxon_lumped )
all.sub.m <- left_join( all.sub, muse )
sub.full <- right_join( sub, all.sub.m )
sub.full$Abundance[ is.na(sub.full$Abundance) ] <- 0
sub.full$Zone <- factor( sub.full$Zone, levels=c('HIGH','MID','LOW'))

ggplot( sub.full, aes(x=Year, y=Abundance, col=Zone)) +  facet_grid(Site~taxon_lumped) +
  geom_point() + geom_smooth() + ggtitle("Non-algal % Cover") +  theme(axis.text.x=element_text(angle=45, hjust=1))

ggplot( sub.full, aes(x=Year, y=Abundance, col=Shore_height_cm, group=Zone)) +  facet_grid(Site~taxon_lumped) +
  geom_point() + geom_smooth() + ggtitle("Non-algal % Cover") +  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_color_viridis()



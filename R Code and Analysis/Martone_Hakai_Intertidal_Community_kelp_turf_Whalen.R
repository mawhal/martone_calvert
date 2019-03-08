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
# Start with West Beach
mwest <- am[am$Site=="West Beach" & am$Zone=="LOW",]
mhigh <- am[ am$Zone!="HIGH",]
muse  <- am
drestrict <- d.simple[ d.simple$UID %in% muse$UID , ] 
dm <- left_join( drestrict, am )



# look at kelp cover across sites and zones
kelp <- dm %>%
  group_by( UID, Date, Quadrat, Meter.point, Site, Zone, Year, Shore_height_cm ) %>%
  filter( Kelp.Fucoid.Turf == "Kelp" ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T) )
# get back abundances for other UIDs in west and fill with zero
kelp.full <- right_join( kelp, muse )
kelp.full$Abundance[ is.na(kelp.full$Abundance) ] <- 0

ggplot( kelp.full, aes(x=Year, y=Abundance)) +  facet_grid(Site~Zone) +
  geom_point() + geom_smooth() + ggtitle("Total Kelp % Cover")



# turfs
turf <- dm %>%
  group_by( UID, Date, Quadrat, Meter.point, Site, Zone, Year, Shore_height_cm ) %>%
  filter( Kelp.Fucoid.Turf %in% c('red turf','brown turf','filament_turf','green turf') ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T) )
# get back abundances for other UIDs in west and fill with zero
turf.full <- right_join( turf, muse )
turf.full$Abundance[ is.na(turf.full$Abundance) ] <- 0

ggplot( turf.full, aes(x=Year, y=Abundance)) +  facet_grid(Site~Zone) +
  geom_point() + geom_smooth() + ggtitle("Total Turf % Cover")

# merge total kelp with total turf
kelp.full$taxon_lumped <- "TOTAL.KELP"
turf.full$taxon_lumped <- "TOTAL.TURF"
kelp.comb <- full_join( kelp.full, turf.full )
kelp.comb$Zone <- factor( kelp.comb$Zone, levels=c('HIGH','MID','LOW' ))

ggplot( kelp.comb, aes(x=Year, y=Abundance, col=taxon_lumped)) +  facet_grid(Site~Zone) +
  geom_point() + geom_smooth() + ggtitle("Kelp + turf trajectories by Site and Zone") +  theme(axis.text.x=element_text(angle=45, hjust=1))

# spread out the data so we can plot kelp against turf
kt.spread <- kelp.comb %>% spread( taxon_lumped, Abundance )

ggplot( kt.spread, aes(x=TOTAL.KELP, y=TOTAL.TURF)) +  facet_grid(Site~Zone, scales="free") +
  geom_point() + geom_smooth(method='lm') + ggtitle("Turf vs Kelp by Site and Zone") 

# remove high plots
kts.high <- kt.spread[ kt.spread$Zone != "HIGH",]
ggplot( kts.high, aes(x=TOTAL.KELP, y=TOTAL.TURF)) +  #facet_grid(Site~Zone, scales="free") +
  geom_point() + geom_smooth(se=F) + geom_smooth(method=lm)+ ggtitle("Turf vs Kelp in MID and LOW") 


# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 26 February 2018

# This script runs ordination procedures on data

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(vegan)
library(psych)



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
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  # filter( non.alga.flag=="Algae" ) %>%
  filter( motile_sessile=="sessile" ) %>%
  group_by( UID, taxon_lumped ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))


# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.simple %>%
  spread( taxon_lumped, Abundance, fill=0 )

# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
muse  <- am
# restrict to rows selected in metadata
d.comm <- d.comm[ d.comm$UID %in% muse$UID, ] 

## ---- non-use if more than algae included
# # there is a quadrat without any metadata, remove this from the metadata
# noalgae <- anti_join( muse, d.comm )
# noalgae$UID
# mclean <- muse[ muse$UID != noalgae$UID, ]



# Sort metadata and community matrix to be the same order
d.comm.order <- d.comm[ order(match(d.comm$UID, muse$UID)),]
cbind( d.comm.order$UID, muse$UID )

# remove UID column from community data
comm <- as.matrix(d.comm.order[,-1])


## ordinate everything separately, and use metadata to label things
metaMDS( comm, distance = "bray", k = 3 )


## ordinate average cover per transect
# merge data with metadata
duse <- d[ d$UID %in% muse$UID, ]
dm <- left_join( duse, muse )
dmean <- dm %>% 
  group_by( Year, Site, Zone, taxon_lumped ) %>%
  summarise( Abundance=mean(Abundance) )
  
# spread out
d.comm.mean <- dmean %>%
  spread( taxon_lumped, Abundance, fill=0 )


# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
comm <- as.matrix(d.comm.mean[,-c(1:3)])

## ordinate means
mds2 <- metaMDS( comm, distance = "bray", k = 2 )


## plotting
par(mar=c(5,4,2,2)+0.1)
plot(mds2,"sites",type="n")
# good segregration along zones
meta$Zone2 <- factor( meta$Zone, levels=c('HIGH','MID','LOW') )
ordispider( mds2,meta$Zone2, col=c('#33a02c','#a6cee3','#1f78b4') )
# ordihull( mds2,meta$Zone, col=c('#33a02c','#a6cee3','#1f78b4') )
ordiellipse( mds2,meta$Zone2, col=c('#33a02c','#a6cee3','#1f78b4') )
points(mds2,pch=21,bg="white")


## look at relationship with elevation
elevm <- dm %>% 
  group_by( Year, Site, Zone ) %>%
  summarise( Shore_height_cm=mean(Shore_height_cm, na.rm=T) )
meta <- left_join( meta, elevm )

ordisurf( mds2, meta$Shore_height_cm,  
          kind = "se", conf = 0.95, label = T )


# sites
par(mar=c(5,4,2,2)+0.1)
plot(mds2,"sites",type="n")
# good segregration along zones
meta$Site2 <- factor( meta$Site )
ordispider( mds2,meta$Site2, col=c('red','blue','black') )
# ordihull( mds2,meta$Zone, col=c('#33a02c','#a6cee3','#1f78b4') )
ordiellipse( mds2,meta$Site2, col=c('red','blue','black') )

points(mds2,pch=21,bg="white")

ordisurf( mds2, meta$Year,  
          kind = "se", conf = 0.95, label = T )

# combine site and year
meta$SiteZone <- with(meta, paste(Site, Zone))
meta$Year2 <- factor(meta$Year)
par(mar=c(5,4,2,2)+0.1)
plot(mds2,"sites",type="n")
zone.col <- c('#33a02c','#a6cee3','#1f78b4')
zone.pch <- c(21,22,25)
site.col <- c('red','blue','black')
year.bg  <- c(NA,NA,NA,NA,NA,NA,'slateblue')
ordisegments( mds2, meta$SiteZone, col=site.col[meta$Site2] )
points( mds2, pch=zone.pch[meta$Zone2], col=site.col[meta$Site2], bg=year.bg[meta$Year2] )


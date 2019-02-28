# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen (adapted in part from Sam Starko's diversity script)
# created 26 Feb 2019


# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(vegan)
library(psych)
library(plotrix)


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read.csv( "Data/R Code/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
# am <- am[ am$Site != "Meay Channel", ]

# remove taxa that are not coutned towards subtratum cover (i.e. mobile invertebrates)
# make it easier by replacing NA values for substratum
ds <- ad
ds$motile_sessile[ is.na(ds$motile_sessile) ] <- "Substratum"
ds <- ds[ ds$motile_sessile!="motile", ]


# remove bare space? Not yet

d <- ds



# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which wil reduce the size of the dataset a bit
# restrict this to seaweeds
d.simple <- d %>%
  filter( non.alga.flag=="Algae" ) %>%
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


# there is a quadrat without any metadata, remove this from the metadata
noalgae <- anti_join( muse, d.comm )
noalgae$UID
mclean <- muse[ muse$UID != noalgae$UID, ]

##Sort metadata and community matrix to be the same order
d.comm.order <- d.comm[ order(match(d.comm$UID, mclean$UID)),]
cbind( d.comm.order$UID, mclean$UID )

# remove UID column from community data
comm <- as.matrix(d.comm.order[,-1])

## Steps
# for each quadrat, calculate richness, Shannon diversity, Simpson Diversity, and ENSPIE
## Total abundance - also used for ENSPIE below
mclean$total.algae <- rowSums( comm )
## shannon
mclean$shannon <- diversity( comm, "shannon" )
## simpson
mclean$simpson <- diversity( comm, "simpson" )
## ENSPIE
# the function
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
} 
prop <- comm/mclean$total.algae
mclean$enspie <- apply( prop, 1, ENSPIE )
# richness
pa <- ifelse( comm>0, 1, 0)
mclean$richness <- rowSums( pa )

# splom for all quadrat summaries
pairs.panels( mclean %>% select(total.algae,richness,shannon,simpson,enspie), scale=T, ellipses = FALSE )
# for each quadrat, calculate total cover of algae -- then use this to calculate Coefficient of Variation


# look at patterns over time
ggplot( mclean, aes(y=richness,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mclean, aes(y=shannon,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mclean, aes(y=simpson,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mclean, aes(y=enspie,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()




# consider the range of variation in estimates over time
divvar <- mclean %>%
  group_by( Year, Site, Zone ) %>%
  summarise( meana=mean(total.algae), ea=sd(total.algae),
             meanr=mean(richness), er=sd(richness),
             meand=mean(shannon), ed=sd(shannon),
             means=mean(simpson), es=sd(simpson),
             meane=mean(enspie), ee=sd(enspie)    ) %>%
  mutate( cva = ea/meana, cvr = er/meanr, cvd = ed/meand, cvs = es/means, cve = ee/meane )

pairs.panels( divvar %>% ungroup() %>% select(cva,cvr,cvd,cvs,cve), scale=T, ellipses = FALSE )
pairs.panels( divvar %>% ungroup() %>% select(meanr,meand,means,meane,cva), scale=T, ellipses = FALSE )

ggplot( divvar, aes(x=meanr,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean transect species richness")

# summarize over time
divvar2 <- mclean %>%
  group_by( Site, Zone ) %>%
  summarise( meana=mean(total.algae), ea=sd(total.algae),
             meanr=mean(richness), er=sd(richness),
             meand=mean(shannon), ed=sd(shannon),
             means=mean(simpson), es=sd(simpson),
             meane=mean(enspie), ee=sd(enspie)    ) %>%
  mutate( cva = ea/meana, cvr = er/meanr, cvd = ed/meand, cvs = es/means, cve = ee/meane )

pairs.panels( divvar2 %>% ungroup() %>% select(meanr,meand,means,meane,cva), scale=T, ellipses = FALSE )

a <- ggplot( divvar2, aes(x=meanr,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean species richness") + geom_smooth() + ylim(c(-0.1,1.2))
b <- ggplot( divvar2, aes(x=meand,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean Shannon diversity") + geom_smooth() + ylim(c(-0.1,1.2))
c <- ggplot( divvar2, aes(x=means,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean Simpson diversity") + geom_smooth() + ylim(c(-0.1,1.2))
d <- ggplot( divvar2, aes(x=meane,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean effective # species") + geom_smooth() + ylim(c(-0.1,1.2))


library( cowplot )
plot_grid( a,b,c,d, ncol=4 )

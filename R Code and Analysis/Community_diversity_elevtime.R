# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen (adapted in part from Sam Starko's diversity script)
# created 26 Feb 2019


# load libraries
library(tidyverse)
library(vegan)
library(psych)
library(plotrix)


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]

# remove taxa that are not coutned towards subtratum cover (i.e. mobile invertebrates)
# make it easier by replacing NA values for substratum
ds <- ad
ds$motile_sessile[ is.na(ds$motile_sessile) ] <- "Substratum"
ds <- ds[ ds$motile_sessile=="sessile", ]
# ds <- ds[ ds$non.alga.flag=="Algae", ]


# remove bare space? Not yet

d <- ds
# # split UID so we can average by transect
splits <- strsplit( as.character(d$UID), " " )
d$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))


# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which wil reduce the size of the dataset a bit
# restrict this to sessile taxa
d.simple <- d %>%
  # filter( motile_sessile=="sessile" ) %>%
  group_by( UID, transect, taxon_lumped ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) %>% 
  arrange(UID)
# # calculate average abundance by transect
# d.trans <- d.simple %>%
#   group_by( transect,taxon_lumped ) %>%
#   summarize( Abundance=mean(Abundance, na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.simple %>%
  spread( taxon_lumped, Abundance, fill=0 )


# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
muse  <- am
splits <- strsplit( as.character(muse$UID), " " )
muse$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))
muse <- arrange(muse,UID)

# restrict to rows selected in metadata
d.comm <- d.comm[ d.comm$transect %in% muse$transect, ] 


# # there is a quadrat without any metadata, remove this from the metadata
# noalgae <- anti_join( muse, d.comm )
# noalgae$UID
# mclean <- muse[ muse$UID != noalgae$UID, ]
mclean <- muse
# define levels for zones
mclean$Zone <- factor( mclean$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
mclean$Site <- factor( mclean$Site, levels = c("West Beach", "Fifth Beach", "North Beach" ))
# define Year factor
# mclean$Year <- factor( mclean$Year, ordered= TRUE )
mclean$transect <- with( mclean, paste(Site,Zone,Year,sep = " ") )
mtrans <- mclean %>%
  group_by( transect, Site, Zone, Year ) %>%
  summarize( Shore_height_cm=mean(Shore_height_cm,na.rm=T) )

##Sort metadata and community matrix to be the same order
# d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
# cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
comm <- as.matrix(d.comm[,-c(1,2)])

anti_join(d.comm[,1:2], mclean[,c("UID", "transect")] )
anti_join(mclean[,c("UID", "transect")] , d.comm[,1:2]  )
## Steps
# for each quadrat, calculate richness, Shannon diversity, Simpson Diversity, and ENSPIE
## Total abundance - also used for ENSPIE below
mclean$total.cover <- rowSums( comm )
## shannon
mclean$shannon <- diversity( comm, "shannon" )
## simpson
mclean$simpson <- diversity( comm, "simpson" )
## ENSPIE
# the function
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
} 
prop <- comm/mclean$total.cover
mclean$enspie <- apply( prop, 1, ENSPIE )
# richness
pa <- ifelse( comm>0, 1, 0)
mclean$richness <- rowSums( pa )

# splom for all quadrat summaries
psych::pairs.panels( mclean %>% ungroup() %>% select(total.cover,richness,shannon,simpson,enspie), 
              scale=F, ellipses = FALSE )
# for each quadrat, calculate total cover of algae -- then use this to calculate Coefficient of Variation


# look at patterns over time
ggplot( mclean, aes(y=simpson,x=Shore_height_cm,col=Year)) + 
  geom_smooth(aes(group=Year), method='loess',se=T,alpha=0.15) + viridis::scale_color_viridis(option = "E")
ggplot( mclean, aes(y=simpson,x=Shore_height_cm,col=Year)) + 
  geom_point(alpha=0.5) + viridis::scale_color_viridis(option = "E")
mean_zone <- mclean %>% group_by(Year,Zone) %>% 
  summarise( enspie=mean(enspie),richness=mean(richness),shannon=mean(shannon), simpson=mean(simpson),Shore_height_cm=mean(Shore_height_cm,na.rm=T) )
ggplot( mclean, aes(y=enspie,x=Shore_height_cm,col=Zone)) + facet_wrap(~Year) +
  geom_point(alpha=0.67) + #col='slateblue',
    # geom_smooth(aes(group=Zone),method='loess',se=T,col='black',method.args = list(family = "gaussian")) +
  geom_smooth(aes(group=1),method='loess',se=T,col='black',method.args = list(family = "gaussian")) +
  geom_point(aes(fill=Zone), data=mean_zone,size=3,shape=21, col='black') +
  xlab("Shore height (m above MLLW)") + ylab("Effective number of species") +
  scale_x_continuous(breaks=c(100,200,300),labels=c(1,2,3)) + 
  scale_color_viridis(discrete=TRUE, begin = 0.9,end = 0) +
  scale_fill_viridis(discrete=TRUE, begin = 0.9,end = 0) +
  # scale_y_continuous(trans='log2') +
  theme_bw()
ggplot( filter(mclean,Year%in%c(2012,2019)), aes(y=enspie,x=Shore_height_cm,group=Year,col=Year)) + 
  geom_point( alpha=0.75 ) +
  geom_smooth(method='loess',se=T,method.args = list(family = "symmetric")) +
    xlab("Shore height (m above MLLW)") + ylab("ENSPIE") +
  scale_x_continuous(breaks=c(100,200,300),labels=c(1,2,3)) + 
  scale_y_continuous(breaks=seq(0,12,2)) +
  theme_bw()

# 
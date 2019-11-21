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
ad <- read.csv( "Data/R code for Data Prep//Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read.csv( "Data/R code for Data Prep//Output from R/Martone_Hakai_metadata.csv" )

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
# # split UID so we can average by transect
splits <- strsplit( as.character(d$UID), " " )
d$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))


# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which wil reduce the size of the dataset a bit
# restrict this to sessile taxa
d.simple <- d %>%
  filter( motile_sessile=="sessile" ) %>%
  group_by( UID, transect, taxon_lumped ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))
# calculate average abundance by transect
d.trans <- d.simple %>%
  group_by( transect,taxon_lumped ) %>%
  summarize( Abundance=mean(Abundance, na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.trans %>%
  spread( taxon_lumped, Abundance, fill=0 )


# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
muse  <- am
splits <- strsplit( as.character(muse$UID), " " )
muse$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))

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
d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
comm <- as.matrix(d.comm.order[,-1])

## Steps
# for each quadrat, calculate richness, Shannon diversity, Simpson Diversity, and ENSPIE
## Total abundance - also used for ENSPIE below
mtrans$total.cover <- rowSums( comm )
## shannon
mtrans$shannon <- diversity( comm, "shannon" )
## simpson
mtrans$simpson <- diversity( comm, "simpson" )
## ENSPIE
# the function
ENSPIE <- function(prop){
  ifelse( sum(prop,na.rm=T)>0, 1 / sum(prop^2, na.rm=T), NA ) 
} 
prop <- comm/mtrans$total.cover
mtrans$enspie <- apply( prop, 1, ENSPIE )
# richness
pa <- ifelse( comm>0, 1, 0)
mtrans$richness <- rowSums( pa )

# splom for all quadrat summaries
pairs.panels( mtrans %>% ungroup() %>% select(total.cover,richness,shannon,simpson,enspie), 
              scale=F, ellipses = FALSE )
# for each quadrat, calculate total cover of algae -- then use this to calculate Coefficient of Variation


# look at patterns over time
ggplot( mtrans, aes(y=richness,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mtrans, aes(y=shannon,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mtrans, aes(y=simpson,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
ggplot( mtrans, aes(y=enspie,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()


ggplot( mtrans, aes(y=enspie,x=Year)) + facet_grid(Site~Zone) + 
  geom_line() + geom_point() + #geom_smooth(se=F) + 
  ylab("Effective number of species")
ggplot( mtrans, aes(y=richness,x=Year)) + facet_grid(Site~Zone) + 
  geom_line() + geom_point() +
  ylab("Species richness")


# make mtrans longer and include both richness and ENSPIE in the same figure
mlong <- mtrans %>%
  select( transect, Site, Zone, Year, Shore_height_cm, enspie, richness ) %>%
  group_by( transect, Site, Zone, Year, Shore_height_cm ) %>%
  gather( Measure, species, -transect, -Site, -Zone, -Year, -Shore_height_cm )

mlong$Measure[mlong$Measure=="enspie"] <- "Effective"
mlong$Measure[mlong$Measure=="richness"] <- "Total"

windows(5,4)
ggplot( mlong, aes(y=species,x=Year, col=Measure)) + facet_grid(Site~Zone) + 
  geom_line() + geom_point() +
  ylab("Number of species") +
  scale_color_manual( values=c("grey","black")) +
  theme_bw() 

# consider the range of variation in estimates over time
divvar <- mtrans %>%
  group_by( Year, Site, Zone ) %>%
  summarise( meana=mean(total.cover), ea=sd(total.cover),
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
divvar2 <- mtrans %>%
  group_by( Site, Zone ) %>%
  summarise( meana=mean(total.cover), ea=sd(total.cover),
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

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


# # look at patterns over time
# ggplot( mclean, aes(y=richness,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=shannon,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=simpson,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()
# ggplot( mclean, aes(y=enspie,x=factor(Year))) + facet_grid(Site~Zone) + geom_boxplot()


ggplot( mclean, aes(y=enspie,x=Year)) + facet_grid(Site~Zone) + 
  geom_line() + geom_point() + #geom_smooth(se=F) + 
  ylab("Effective number of species")
ggplot( mclean, aes(y=richness,x=Year)) + facet_grid(Site~Zone) + 
  geom_line() + geom_point() +
  ylab("Species richness")


# make mclean longer and include both richness and ENSPIE in the same figure
mlong <- mclean %>%
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
divvar <- mclean %>%
  group_by( Site, Zone ) %>%
  summarise( meana=mean(total.cover), ea=sd(total.cover),
             meanr=mean(richness), er=sd(richness),
             meand=mean(shannon), ed=sd(shannon),
             means=mean(simpson), es=sd(simpson),
             meane=mean(enspie), ee=sd(enspie),
             elev = mean(Shore_height_cm,na.rm=T) ) %>%
  mutate( cva = ea/meana, cvr = er/meanr, cvd = ed/meand, cvs = es/means, cve = ee/meane )

ggplot( divvar, aes(x=meanr,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean transect species richness")
ggplot( divvar, aes(x=elev,y=cva,size=meanr,col=Site) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean transect elevation (cm)")

# collapse and plot all together
divplot <- divvar %>% 
  group_by(Site, Zone, cva) %>% 
  gather(key="estimate",value="mean", meanr, meand, means, meane )

ggplot( data=divplot, aes(x=mean, y=cva)) + facet_wrap(~estimate, scales="free") + 
  geom_smooth(method='lm',se=F) +
  geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(col=Site,size=Zone))



## above calculations include variation across and within times, exclude this variation within time?
yearly <- mclean %>% 
  group_by( Site, Zone, Year ) %>% 
  summarize( meana=mean(total.cover), 
             meanr=mean(richness),
             meand=mean(shannon),
             means=mean(simpson),
             meane=mean(enspie) )

# summarize over time
divvar2 <- yearly %>%
  group_by( Site, Zone ) %>%
  summarise( gmeana=mean(meana), ea=sd(meana),
             gmeanr=mean(meanr), er=sd(meanr),
             gmeand=mean(meand), ed=sd(meand),
             gmeans=mean(means), es=sd(means),
             gmeane=mean(meane), ee=sd(meane)    ) %>%
  mutate( cva = ea/gmeana, cvr = er/gmeanr, cvd = ed/gmeand, cvs = es/gmeans, cve = ee/gmeane )

# collapse and plot all together
divplot2 <- divvar2 %>% 
  group_by(Site, Zone, cva) %>% 
  gather(key="estimate",value="mean", gmeanr, gmeand, gmeans, gmeane )

ggplot( data=divplot2, aes(x=mean, y=1/cva)) + facet_wrap(~estimate, scales="free") + 
  geom_smooth(method='lm',se=T) +
  geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(col=Site,size=Zone))

summary(lm(mean~cva, data=filter(divplot2,estimate=="gmeanr")))
ggplot( data=filter(divplot2,estimate=="gmeanr"), aes(x=mean, y=1/cva)) + 
  geom_smooth(method='lm',se=T) +
  geom_smooth(aes(group=Zone,lty=Zone),method='lm',se=F, col='black') +
  geom_point(aes(fill=Site,size=Zone),pch=21) +
  ylab( expression(paste("Cover stability (",mu,"/",sigma,")")) ) + xlab("Mean species richness") +
  viridis::scale_fill_viridis(discrete=T) +
  theme_classic()




a <- ggplot( divvar2, aes(x=gmeanr,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean species richness") + geom_smooth() + ylim(c(-0.1,1.2))
b <- ggplot( divvar2, aes(x=gmeand,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean Shannon diversity") + geom_smooth() + ylim(c(-0.1,1.2))
c <- ggplot( divvar2, aes(x=gmeans,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean Simpson diversity") + geom_smooth() + ylim(c(-0.1,1.2))
d <- ggplot( divvar2, aes(x=gmeane,y=cva) ) + geom_point() +
  ylab( "CV( total algal % cover )") + xlab("Mean effective # species") + geom_smooth() + ylim(c(-0.1,1.2))


library( cowplot )
plot_grid( a,b,c,d, ncol=4 )

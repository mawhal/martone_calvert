# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# started 13 May 2019

# This script uses looks at correlations between temperature and community attributes (including biodiversity metrics)

# load libraries
library( tidyverse )
library( vegan )


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_lump_function.csv", stringsAsFactors = FALSE )
# all metadata
am <- read.csv( "Data/R Code/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )
# temperature data
sst <- read_csv( "Data/Environmental Data/PineIsland_summary.csv" )




## Data cleaning for Analysis -- consider moving part of this to another script
# # remove 2011 data
muse <- am[ am$Year != "2011", ]
# # remove Meay Channel
# ## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
# ## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
# muse <- muse[ muse$Site != "Meay Channel", ]
# muse <- droplevels(muse)
# muse <- am

# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to algae only
d <- dm %>% 
  filter( non.alga.flag =="Algae" )
bare <- dm %>%
  filter( taxon_lumped2 == "Bare rock" )

# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  group_by( UID, Year, Site, Zone, taxon_lumped2 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))

# average cover per transect
dmean <- d.simple %>% 
  group_by( Year, Site, Zone, taxon_lumped2 ) %>%
  summarise( Abundance=mean(Abundance) )
# average bare rock cover
bmean <- bare %>%
  group_by( Year, Site, Zone ) %>%
  summarize( bare = mean(Abundance))

# get average total algal cover
dtotal <- d.simple %>%
  group_by( Year, Site, Zone, UID ) %>%
  summarise( total.cover = sum(Abundance) ) %>%
  group_by( Year, Site, Zone ) %>%
  summarise( mean.cover = mean(total.cover), sd.cover=sd(total.cover) )


# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- dmean %>%
  spread( taxon_lumped2, Abundance, fill=0 )



# isolate the community, site, and sample data
comm <- d.comm[,-c(1:3)]
# site <- d.comm$Site
# tran <- apply( d.comm[,c(2,3)],1,paste, collapse="." )
# year <- d.comm$Year
# year2 <- as.character(year)
# year2[year2!=2016] <- ""
# yearmod <- factor(year,ordered=T)
# zone <- d.comm$Zone

# calculate alpha diversity scores -- revisit the possible need to think of the transect as the local level of diversity rather than quadrat
alpha <- diversity( comm, index = "shannon" )


# combine total cover with shannon diversity
d.sum <- data.frame( dtotal, alpha )
d.sum <- left_join( d.sum, bmean )

# add sst
d.sst <- left_join( d.sum, sst )





# make a few figures
d.sst$zoneyear <- with(d.sst, paste(Zone,Year,sep="-") )


ggplot( d.sst, aes(x=Anomaly,y=bare,col=Year)) + 
  facet_wrap( ~Site ) +
  geom_point(size=3) +
  geom_smooth(aes(lty=Zone),method='lm', se=F, col='black')

ggplot( d.sst, aes(x=Year,y=bare,col=Year)) + 
  facet_wrap( ~Site ) +
  geom_point(size=3) +
  geom_smooth(aes(lty=Zone), se=F, col='black')

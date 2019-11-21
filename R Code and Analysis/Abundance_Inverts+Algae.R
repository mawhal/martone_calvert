# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# created 25 September 2019 -- adpated from Species_TimeSeries

# This script produces figures of the density of algae and invertebrates, saving figures as pdf


# load libraries
library(tidyverse)
library(viridis)

## read data files
# all data
ad <- read.csv( "Data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv")
# all metadata
am <- read.csv("Data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )




# taxa to remove...for now
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]



## Customizations to carry through to figures

# Choose the group
dfilt <- ad %>%
  filter( motile_sessile == "sessile" ) 
  

# merge data and metadata
# get all instances of a particular taxon
# to include all quadrats us full_join, or use left_join for quads with the taxon
d <- full_join( dfilt, am )
d$Abundance[ is.na( d$Abundance) ] <- 0

# remove Meay Channel, and 2011
d <- d %>%
  filter( Site != "Meay Channel", Year != 2011 )


# make abundances numeric
sort(unique(d$Abundance))
d$Abundance <- as.numeric( d$Abundance )
# define leveles for zones
d$Zone <- factor( d$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
d$Site <- factor( d$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))
# define Year factor
d$Year <- factor( d$Year, ordered= TRUE )


# because we are focused on comparing two groups rather that the individual taxa, we will sum by group
d <- d %>%
  group_by( UID, Site, Zone, Year, Quadrat, non.alga.flag ) %>%
  summarize( Abundance=sum(Abundance) )

# set West Beach as the first level
d$Site <- relevel(d$Site, ref="West Beach") 
# rename animal to sessile invertebrate
d$non.alga.flag <- as.character(d$non.alga.flag)
d$non.alga.flag[d$non.alga.flag=="Animal"] <- "Sessile\nInvertebrate"

# all sites
# time trends in different tidal zones
windows(5,6)
(ggzone <- ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance, col=non.alga.flag)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=F ) +
    stat_summary( fun.data = "mean_cl_boot", size = 0.5 ) +
    stat_summary( fun.y = "mean", geom="line", size = 0.5 ) +
    # geom_point( alpha=0.4,col='slateblue' ) + ggtitle( taxon ) + 
    xlab("Year") +
    scale_x_continuous(breaks = seq(2010,2018,2) ) ) +
    scale_color_manual( values=c("slateblue4","slateblue1") ) +
    theme_bw() + theme(legend.title=element_blank())
ggsave( "R Code and Analysis/Figs/algae+invert_means.pdf", ggzone, "pdf" )

# make a stacked version
# calculate means
windows(5,4)
dmean <- d %>%
  group_by( Site, Zone, Year, non.alga.flag ) %>%
  summarize( Abundance = mean(Abundance,na.rm=T) )
(ggstack <- ggplot( dmean, aes(x=as.numeric(as.character(Year)),y=Abundance, fill=non.alga.flag)) + 
    facet_grid(Site~Zone) + 
    # geom_smooth( se=F ) +
    geom_area( ) + 
    xlab("Year") +
    scale_x_continuous(breaks = seq(2010,2018,2) ) ) +
    scale_fill_manual( values=c("slateblue1","slateblue4") ) +
  theme_bw() + theme(legend.title=element_blank()) 


# subset of sites where elevation has been measured
delev <- d[ d$Site != "Meay Channel", ]
windows(10,4)
(ggheight <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance)) + 
    facet_grid(Site~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="quasipoisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) )
# ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation_wide.pdf"), ggheight, "pdf" )

windows(4,6)
(ggheight2 <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance,group=Year,col=Year )) + 
    facet_wrap(~Site,ncol=1, scales = "free_y") + 
    geom_point(alpha=0.75) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="poisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) +
    # geom_smooth(aes(group=1)) + 
    scale_x_continuous(trans='log10') ) +
  scale_color_viridis_d( direction=-1 )

# ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation.pdf"), ggheight2, "pdf" )


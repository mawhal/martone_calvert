# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script produces figures of the density of a chosen taxa, saving figures as pdf
taxon <- "Alaria"
# sampler <- "Sandra" --- figure out how to add a switch here that we can add to filenames

# set options
options(stringsAsFactors = FALSE)

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

# Choose a taxon 
sort( unique( ad$taxon_lumped ))
# use general exp to pull several groups if needed
sort(unique( ad$taxon_lumped[ grep( paste0(taxon,"*"), ad$taxon_lumped ) ]  ))
dtax   <- ad[ grep( paste0(taxon,"*"), ad$taxon_lumped ), ]
#dtax  <-  dtax[ -grep( "Ectocarpus*",dtax$Taxon), ]


# merge data and metadata
# get all instances of a particular taxon
# to include all quadrats us full_join, or use left_join for quads with the taxon
d <- full_join( dtax, am )
d$Abundance[ is.na( d$Abundance) ] <- 0





# *** to add *** include full names for each unique instance of Sampler


# # chose an observer/recorder
# sort(unique( am$Sampler ))
# sort(unique( am$Recorder ))
# data.frame( Sampler=sort(unique( c(am$Sampler,am$Recorder))), fullName="_________"  )
# # Sandra C Lindstrom
# sandra     <- c( "Sandra", "SCL", "SL" )
# # select the plots we want from the metadata
# dperson <- d %>%
#   filter( Sampler %in% sandra | Recorder %in% sandra )
# d<- dperson


# make abundances numeric
sort(unique(d$Abundance))
d$Abundance <- as.numeric( d$Abundance )
# define leveles for zones
d$Zone <- factor( d$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
d$Site <- factor( d$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))
# define Year factor
d$Year <- factor( d$Year, ordered= TRUE )

# all sites
# time trends in different tidal zones
windows(5,6)
(ggzone <- ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=TRUE, col='black' ) +
    stat_summary( fun.data = "mean_cl_boot", colour = "slateblue4", size = 0.5 ) +
    stat_summary( fun.y = "mean", geom="line", colour = "slateblue4", size = 0.5 ) +
    geom_point( alpha=0.4,col='slateblue' ) + ggtitle( taxon ) + 
    xlab("Year") )

ggsave( paste0("R Code and Analysis/Figs/",taxon,"_zone.pdf"), ggzone, "pdf" )

# subset of sites where elevation has been measured
delev <- d[ d$Site != "Meay Channel", ]
windows(10,4)
(ggheight <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance)) + 
    facet_grid(Site~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="quasipoisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) )
ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation_wide.pdf"), ggheight, "pdf" )

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

ggsave( paste0("R Code and Analysis/Figs/",taxon,"_elevation.pdf"), ggheight2, "pdf" )

# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 29 January 2019

# This script produces figures of the density of a chosen taxa, saving figures as pdf
taxon <- "Fucus"
# sampler <- "Sandra" --- figure out how to add a switch here that we can add to filenames

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(viridis)

## read data files
# all data
ad <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data.csv")
# all metadata
am <- read.csv("Data/R Code/Output from R/Martone_Hakai_metadata.csv" )


## Deal with trace cover and other oddities
# # replace all commas with periods for Abundance
# ad$Abundance <- as.numeric( gsub( ",", "[.]", ad$Abundance ) )
# change "present" to 0,5
ad$Abundance <- gsub( "present", "0.5", ad$Abundance )
# change 1 of Fucus to "trace"
ad$Abundance <- gsub( "1 on Fucus", "trace", ad$Abundance ) 
# change trace to 0.5% cover
ad$Abundance <- gsub( "t.*", "0.5", ad$Abundance, ignore.case = TRUE ) 

# accept only the first thing if separated by certain characters
# asplit <- strsplit( ad$Abundance, split =c("/|;"))
# ad$Abundance <- unlist( lapply( asplit, function(z) z[1] ))
sort(unique(ad$Abundance))
# things to fix

# spaces
# periods after numbers -- e.g. "0.5."

# taxa to remove...for now
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]



## Customizations to carry through to figures

# Choose a taxon 
sort( unique( d$Taxon ))
# use general exp to pull several groups if needed
sort(unique( d$Taxon[ grep( paste0(taxon,"*"), d$Taxon ) ]  ))
dtax   <- ad[ grep( paste0(taxon,"*"), ad$Taxon ), ]
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
windows(5,7)
(ggzone <- ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    geom_point(alpha=0.2) + ggtitle( taxon ) + 
    geom_smooth( se=TRUE ) +
    xlab("Year") )

ggsave( paste0("Figs/",taxon,"_zone.pdf"), ggzone, "pdf" )

# subset of sites where elevation has been measured
delev <- d[ d$Site != "Meay Channel", ]
windows(10,4)
(ggheight <- ggplot( delev, aes(x=Shore_height_cm,y=Abundance)) + 
    facet_grid(Site~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( taxon ) + 
    geom_smooth(method="glm", method.args=list(family="poisson"), 
                formula = ceiling(y) ~ poly(x,2), 
                se=FALSE, lwd=0.5) )
ggsave( paste0("Figs/",taxon,"_elevation_wide.pdf"), ggheight, "pdf" )

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

ggsave( paste0("Figs/",taxon,"_elevation.pdf"), ggheight2, "pdf" )

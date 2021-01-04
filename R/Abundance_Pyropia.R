# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# created 28 Feb 2019

# This script produces shows patterns of Pyropia abbottiae and other Pyropia species

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(viridis)

## read data files
# all data
ad <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_lump_function.csv")
# all metadata
am <- read.csv("Data/R Code/Output from R/Martone_Hakai_metadata.csv" )




# taxa to remove...for now
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]


## Identitfy customizations to carry through to figures
# Choose a taxon 
sort( unique( ad$Taxon ))
sort( unique( ad$taxon_lumped ))
# use general exp to pull several groups if needed
sort(unique( ad$Taxon[ grep( paste0('Pyropia',"*"), ad$taxon_lumped ) ]  ))
dpyr    <- ad[ grep( paste0('Pyropia',"*"), ad$taxon_lumped ), ]
sort(unique(dpyr$taxon_lumped))


# merge data and metadata
# include all cases?
d <- left_join( dpyr,am )


# write to disk
write.csv( d, "R Code and Analysis/output from r/Pyropia.csv", row.names=FALSE )


# define leveles for zones
d$Zone <- factor( d$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
d$Site <- factor( d$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))

# all sites
# time trends in different tidal heights
windows(5,7)
(ggzone <- ggplot( d, aes(x=Year,y=Abundance)) + facet_grid(Site~Zone, scales="free_y") + 
    geom_point(alpha=0.2) + geom_smooth() )

ggsave( paste0("Figs/AlariaSaccKaty_",site,"_zone.pdf"), ggzone, "pdf" )



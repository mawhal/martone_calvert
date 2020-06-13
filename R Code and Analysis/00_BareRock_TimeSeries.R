# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# library(profvis)
# profvis({
# This script produces figures of the density of a chosen taxa, saving figures as pdf
taxon <- "Bare rock"

# taxa with predicted increases over time
# taxon <- "Mytilus"
# # taxon <- "Lithophyllum"
# # taxon <- "Neopolyporolithon reclinatum"
# # taxon <- "Lithothamnion phymatodeum"
# # taxon <- "Dilsea.californica"
# # sampler <- "Sandra" #--- figure out how to add a switch here that we can add to filenames
# 
# taxon <- "Alaria"

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




# taxa to remove
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]

# customize sites and years
years <- 2012:2019
sites <- c("Fifth Beach", "North Beach", "Foggy Cove")
metause <- am %>% 
  filter( Year %in% years ) %>% 
  filter( Site %in% sites )

## Customizations to carry through to figures

# Choose a taxon 
sort( unique( ad$taxon_lumped2 ))
# use general exp to pull several groups if needed
sort(unique( ad$taxon_lumped2[ grep( paste0(taxon,"*"), ad$taxon_lumped2 ) ]  ))
dtax   <- ad[ grep( paste0(taxon,"*"), ad$taxon_lumped2 ), ]
#dtax  <-  dtax[ -grep( "Ectocarpus*",dtax$Taxon), ]


# merge data and metadata
# get all instances of a particular taxon
# to include all quadrats us full_join, or use left_join for quads with the taxon
d <- left_join( metause, dtax )
d$Abundance[ is.na( d$Abundance) ] <- 0
# # get rid of bad merging
# d <- d %>% 
#   filter( !is.na(Site), !is.na(Zone) )


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
d$Site <- factor( d$Site, levels = sites )
# define Year factor
d$Year <- factor( d$Year, ordered= TRUE )

write_csv(d,"R Code and Analysis/output from r/bare.csv")

# all sites
# time trends in different tidal zones
windows(5,5)
(ggzone <- ggplot( d, aes(x=as.numeric(as.character(Year)),y=Abundance)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    geom_smooth( se=TRUE, col='black' ) +
    geom_point( alpha=0.4,col='slateblue' ) +  
    xlab("Year") + ylab("% primary substrate") +
    ylim( c(0,100) ) +
    scale_x_continuous(breaks = seq(2010,2018,2) ) )

ggsave( paste0("R Code and Analysis/Figs/bare_site.pdf"), ggzone, "pdf" )

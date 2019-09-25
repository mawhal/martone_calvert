# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 28 January 2019

# This script runs diversity analyses on the Rocky Shore Data using package 


# TO DO FOR DATA COMBINE AND CLEAN SCRIPTS
# add together taxa that are not unique to each quadrat
# Change column name for Quadrat.No to Quadrat in the metadata
# Add column to metadata for a unique ID for each quadrat in each zone at each site in each year


# set options
# lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
# options(stringsAsFactors = FALSE)

# load libraries
library(mobr)
# library(dplyr)
# library(tidyverse)

# data(inv_comm) # Community matrix
# data(inv_plot_attr) # Plot attributes data.frame
# inv_mob_in <- make_mob_in(inv_comm, inv_plot_attr)
# inv_mob_in

## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
comm <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_community.csv",stringsAsFactors = FALSE )
# all metadata
am <- read.csv( "Data/R Code/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]
# currently, we do not have geolocation info for North Beach, so remove this one too
am <- am[ am$Site != "North Beach", ]
# try only low zone for now
am <- am[ am$Zone == "LOW", ]

am <- droplevels(am)




### Try to use mobr to show  diversity patterns


# prepare the data
row.names(comm)
row.names(am) <- 1:nrow(am)
sea_mob_in <- make_mob_in( comm , am, c("Long","Lat"), latlong=TRUE )



# exploratory analysis - samples added with increasing spatial proximity
windows()
cols=c( "slateblue",   "darkorange" )  #"black",
plot_rarefaction( sea_mob_in, "Site", "spat", lwd=4, col=scales::alpha(cols,0.8) ) 
# West beachtends to have higher richness at all scales, but especially coarser scales

windows(7,5)
par(mfrow=c(1,2))
plot_rarefaction( sea_mob_in, 'Site', 'indiv', pooled=F, lwd=2,
                 leg_loc='topright', col=scales::alpha(cols,0.1) )

plot_rarefaction( sea_mob_in, 'Site', 'indiv', pooled=T, lwd=4,
                 leg_loc=NA, col=scales::alpha(cols,0.8) )


# Species Abundance Distribution
par(mfrow=c(1,2))
plot_abu( sea_mob_in, 'Site', type='rad', pooled=F, log='x', col=scales::alpha(cols,0.1) )
plot_abu( sea_mob_in, 'Site', type='sad', pooled=F, log='x', col=scales::alpha(cols,0.1) )
plot_abu( sea_mob_in, 'Site', type='rad', pooled=T, log='x', col=scales::alpha(cols,0.8))
# Slightyly higher evenness at West Beach than Fifth Beach

sea_stats <- get_mob_stats( sea_mob_in, group_var = "Site",
                            n_perm = 200 )
windows(9,3.25)
plot(sea_stats, 'S', col=scales::alpha(cols,0.6) )
windows(9,3.25)
plot(sea_stats, 'N', col=scales::alpha(cols,0.6) )
windows(9,3.25)
plot(sea_stats, 'S_n', col=scales::alpha(cols,0.6) )
windows(9,3.25)
plot(sea_stats, 'S_PIE', col=scales::alpha(cols,0.6) )

# West Beach more diverse than Fifth Beach at every scale
#  driven by both richness and abundance
#  although richness difference not significant
#  only significant differences are for abundance, beta of ENSPIE



### Multiscale analysis
sea_deltaS = get_delta_stats( sea_mob_in, 'Site', ref_group='Fifth Beach',
                              type='discrete', log_scale=FALSE, n_perm = 20)
windows( 7,4 )
plot( sea_deltaS, 'West Beach', 'Fifth Beach', display='rarefaction')
plot( sea_deltaS, 'West Beach', 'Fifth Beach', display='delta S')

windows( 7,4 )
plot( sea_deltaS, 'West Beach', 'Fifth Beach', display='ddelta S')

windows( 7,4 )
par(mfrow=c(1,2))
overlap_effects(sea_deltaS, 'West Beach', display='raw', leg_loc = 'bottomright')
overlap_effects(sea_deltaS, 'West Beach', display='stacked', prop=T, leg_loc = NA)

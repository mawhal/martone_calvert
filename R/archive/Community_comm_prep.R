# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen


# This script prepares community data for analysis

# TO DO FOR DATA COMBINE AND CLEAN SCRIPTS
# add together taxa that are not unique to each quadrat
# Change column name for Quadrat.No to Quadrat in the metadata
# Add column to metadata for a unique ID for each quadrat in each zone at each site in each year



# load libraries
# library(mobr)
# library(dplyr)
library(tidyverse)

# data(inv_comm) # Community matrix
# data(inv_plot_attr) # Plot attributes data.frame
# inv_mob_in <- make_mob_in(inv_comm, inv_plot_attr)
# inv_mob_in

## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]
am <- droplevels(am)

# remove taxa that are not counted towards subtratum cover (i.e. mobile invertebrates)
ds <- ad
# ds <- ds[ ds$motile_sessile!="motile", ]

# for now, restrict community analysis to algae only
# ds <- ds %>% 
#   filter( non.alga.flag =="Algae" )

# remove bare space, because we want to look for differences in numbers of individuals
ds <- ds[ ds$Taxon != "bare rock",]  # need to make sure using clean names here

d <- ds

d %>% filter( is.na(taxon_lumped2) ) %>% select(Taxon) %>% unique()

# add up all taxa that appear more than once in a single quadrat (e.g. barnacles or Hildenbrandia) -- go back to datasheets on some of these?
d[ duplicated(d), ] # generalize this to look at particular columns

# spread out all of the community data so sample (site,height,year,quadrat) are rows and species are columns
# add together taxa that are not unique to each quadrat
d.simple <- d %>%
  group_by( UID, taxon_lumped2 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))
d.simple$taxon_lumped2 <- gsub(" ",".",d.simple$taxon_lumped2)


# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.simple %>%
  tidyr::spread( taxon_lumped2, Abundance, fill=0 )


am$UID[ !(am$UID %in% d.comm$UID) ] 
d[ d$UID %in% am$UID[ !(am$UID %in% d.comm$UID) ], ]
# restrict to rows selected in metadata
d.comm <- d.comm[ d.comm$UID %in% am$UID, ] 
# d.comm <- d.comm[,-1]

View(d.comm[,c("UID","Savoiea robusta")])

# write the community data to disk
write_csv( d.comm, "R/output/community.csv" )




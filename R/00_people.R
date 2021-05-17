# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# library(profvis)
# profvis({
# This script produces figures of the density of a chosen taxa, saving figures as pdf

sampler <- "Sandra" #--- figure out how to add a switch here that we can add to filenames

taxon <- "Alaria"

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
sites <- c("Fifth Beach", "North Beach","Foggy Cove", "Meay Channel")
metause <- am %>% 
  filter( Year %in% years ) %>% 
  filter( Site %in% sites )

## Customizations to carry through to figures

# Choose a taxon 
sort( unique( ad$taxon_lumped ))
sort( unique( ad$taxon_revised ))
# use general exp to pull several groups if needed
sort(unique( ad$taxon_revised[ grep( paste0(taxon,"*"), ad$taxon_revised ) ]  ))
dtax   <- ad[ grep( paste0(taxon,"*"), ad$taxon_revised ), ]
#dtax  <-  dtax[ -grep( "Ectocarpus*",dtax$Taxon), ]


# merge data and metadata
# get all instances of a particular taxon, and code missing as zero
# include all quadrats
d <- left_join( metause, dtax )
d$Abundance[ is.na( d$Abundance) ] <- 0
# # get rid of bad merging
# d <- d %>% 
#   filter( !is.na(Site), !is.na(Zone) )
# # query d for consistent plots revisited
# drevisited <- d %>% 
#   group_by( Site, Zone, Meter.point ) %>% 
#   summarize( n = length(Quadrat) ) %>% 
#   arrange( -n )
# dfilt <- left_join(d, drevisited) %>% 
#   filter( n > 5 )
# d <- dfilt

# *** to add *** include full names for each unique instance of Sampler


# chose an observer/recorder
sort(unique( am$Sampler ))
sort(unique( am$Recorder ))
allpeeps <- data.frame( Person = sort(unique( c(am$Sampler,am$Recorder))), fullName=""  )
write_csv(allpeeps, "R/output/all_samplers_recorders.csv")
# combine 
order_of_appearance <- am %>% 
  mutate( year = lubridate::year(Date) ) %>% 
  select(year, Sampler, Recorder) %>% 
  group_by(year) %>% 
  pivot_longer(!year, names_to = "job", values_to = "name" ) %>% 
  select(year, name) %>% 
  distinct() %>% 
  ungroup() %>% group_by(name) %>% 
  summarize( first = min(year), n = length(year) ) %>% 
  arrange(first)

# Sandra C Lindstrom
sandra     <- c( "Sandra", "SCL", "SL" )
# select the plots we want from the metadata
dperson <- d %>%
  filter( Sampler %in% sandra | Recorder %in% sandra )
d<- dperson


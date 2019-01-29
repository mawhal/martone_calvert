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
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(mobr)
data(inv_comm) # Community matrix
data(inv_plot_attr) # Plot attributes data.frame
inv_mob_in <- make_mob_in(inv_comm, inv_plot_attr)
inv_mob_in

## read data files
# all data
ad <- read.csv( "../Data/R Code/Output from R/Martone_Hakai_data.csv")
# all metadata
am <- read.csv("../Data/R Code/Output from R/Martone_Hakai_metadata.csv" )
# add a column for quadrat + SiteHEightYear
am <- am %>%
  mutate( quadSiteHeightYear = paste(Quadrat.No.,SiteHeightYear),
          x=0, y=0 )


# start by filtering out taxa that are not seaweeds


# for now, remove bare space, because we want to look for differences in numbers of individuals
ad <- ad[ ad$Taxon != "bare rock",]  # need to make sure using clean names here



# add up all taxa that appear more than once in a single quadrat (e.g. barnacles or Hildenbrandia) -- go back to datasheets on some of these?
ad[ duplicated(ad), ] # generalize this to look at particular columns

# spread out all of the community data so sample (site,height,year,quadrat) are rows and species are columns
# add together taxa that are not unique to each quadrat
ad.simple <- ad %>%
  mutate( quadSiteHeightYear = paste(Quadrat,SiteHeightYear)) %>%
  select( quadSiteHeightYear, Taxon, Abundance ) %>%
  group_by( quadSiteHeightYear, Taxon ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))

ad.comm <- ad.simple %>%
  spread( Taxon, Abundance, fill=0 )

ad.comm <- as.matrix(ad.comm[,-1])




### Try to use mobr to show some diversity patterns

# prepare the data
sea_mob_in <- make_mob_in( ad.comm , am )

# exploratory analysis - samples added with increasing spatial proximity
windows()
plot_rarefaction( sea_mob_in, "Site", "spat", lwd=4 ) 
# North Beach and Meay Channel have highest richness
# followed by West Beach then Fifth Beach, although West Beach accumulates quickly

par(mfrow=c(1,2))
plot_rarefaction( sea_mob_in, 'Site', 'indiv', pooled=F, lwd=2,
                 leg_loc='topright')
# Meay Channel has an outrageous number of individuals
plot_rarefaction( sea_mob_in, 'Site', 'indiv', pooled=T, lwd=4,
                 leg_loc=NA)
# North Beach again with highest richness, 
# but Meay Channel has much lower richness at given number of individuals
# suggests Meay has lower evenness

# Species Abundance Distribution
par(mfrow=c(1,2))
plot_abu( sea_mob_in, 'Site', type='rad', pooled=F, log='x')
plot_abu( sea_mob_in, 'Site', type='rad', pooled=T, log='x')
# very low evenness at Meay Channel, but this could be due to fewer samples overal

# multiple scales
sea_stats <- get_mob_stats(sea_mob_in, group_var = "Site",
                           n_perm = 200)

sea_stats <- inv_stats

plot(sea_stats, 'S')
plot(sea_stats, 'N')
plot(sea_stats, 'S_n')
plot(sea_stats, 'S_n')
  

# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script looks up tax for different purposes

# load libraries
library( tidyverse )


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv", stringsAsFactors = FALSE )
# all metadata
am <- read.csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )






#### find all taxa from mid shore Fifth Beach #####
fifth_mid_taxa <- ad %>% 
  separate(UID, into = c('Site','Beach','Zone','Year','Quad', sep=" ")) %>% 
  filter( Site ==  "Fifth",  Zone == "MID") %>% 
  filter( non.alga.flag == "Algae" ) %>% 
  select( taxon = taxon_revised, functional_group = new_cat_simple1 ) %>% 
  distinct() %>% 
  arrange( functional_group, taxon )

write_csv( fifth_mid_taxa, "../../../Downloads/Fifth_Beach_MID_taxa.csv")


#### find all taxa from mid and low shore North Beach #####
north_midlow_taxa <- ad %>% 
  separate(UID, into = c('Site','Beach','Zone','Year','Quad', sep=" ")) %>% 
  filter( Site ==  "North",  Zone %in% c("MID", "LOW")) %>% 
  filter( non.alga.flag == "Algae" ) %>% 
  select( taxon = taxon_revised, functional_group = new_cat_simple1 ) %>% 
  distinct() %>% 
  arrange( functional_group, taxon )

write_csv( north_midlow_taxa, "../../../Downloads/North_Beach_MIDLOW_taxa.csv")

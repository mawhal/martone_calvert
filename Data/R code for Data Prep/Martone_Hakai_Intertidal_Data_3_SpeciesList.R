###################################################
### Martone Lab Rocky Intertidal Community Data ###
###     All Data collected on Calvert Island    ###
###################################################

### The purpose of this script is to relate
### species in the raw data 
###   (as it was entered from field sheets)

###################################################
### Data are percent cover of seaweeds from 
###   quadrats along transects in low, mid,
###   and high tide zones
###################################################
###
### code created by Matt Whalen, 16 April 2018
### updated 20 December 2018


## load libaries

library( tidyverse )


## read files

# Data
data <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data.csv" )

#metadata
meta <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

# Unique species from the Data
# sort( unique( data$Taxon ) )
# write this list of unique names to file
write_csv( data.frame(taxon=sort(unique( data$Taxon ))), 
           "data/R Code for Data Prep/Output from R/Martone_Hakai_uniqueTaxa.csv" )



# functional group data
functional <- read_csv("Data/taxa/functional_groups.csv")

# are all taxa representing in both the corrected taxon list and in the functional groups sheet?
corrected_taxa <- read.csv( "data/taxa/CorrectedTaxonList_2019.csv" ) %>% 
  select( -taxon ) %>% distinct()
which( duplicated(corrected_taxa$taxon_corrected))

# merge data with corrected taxa
data <- left_join( data, corrected_taxa, by = c("Taxon" = 'taxon_corrected') )
## merge functional traits with lumped species
data.funct  <- left_join( data, functional, by = c("taxon_lumped3"="taxon"))
# # which lines are messed up?
# extras <- data.funct[ duplicated( data.funct[,c(1:7)] ), ]
# sort( unique(extras$Taxon) )
( no_functional_group <- data.funct %>% filter(motile_sessile=="sessile") %>% filter(is.na(funct_2021)) %>% 
  filter(non.alga.flag == "Algae") %>% 
  select(taxon_lumped3, flag) %>% distinct() )
write_csv( no_functional_group, "Data/taxa/taxa_lacking_functional_grouping.csv" )



# write to disk
write_csv( data.funct, "data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )

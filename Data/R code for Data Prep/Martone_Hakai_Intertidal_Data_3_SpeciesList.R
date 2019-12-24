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
meta<-read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

# Unique species from the Data
sort( unique( data$Taxon ) )
# write this list of unique names to file
write_csv( data.frame(taxon=sort(unique( data$Taxon ))), 
           "data/R Code for Data Prep/Output from R/Martone_Hakai_uniqueTaxa.csv" )

# Load lumping data -- some species are indistinguishable in the field, or were not at the time of the start of the project
lump <- read_csv("data/taxa/TaxonList_corrected_lumped_unique.csv")

# functional group data
functional <- read_csv("Data/taxa/Algae_functional_groups.csv")




## merge data with lumped names
lumped.data <- left_join( data, lump, by = c("Taxon"="taxon_corrected") )
# # which lines are messed up?
# extras <- lumped.data[ duplicated( lumped.data[,c(1:3)] ), ]


## merge functional traits with lumped species
data.funct  <- left_join( lumped.data, functional, by = c("taxon_lumped"="taxon"))
# # which lines are messed up?
# extras <- data.funct[ duplicated( data.funct[,c(1:7)] ), ]
# sort( unique(extras$Taxon) )




# write to disk
write_csv( data.funct, "data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )

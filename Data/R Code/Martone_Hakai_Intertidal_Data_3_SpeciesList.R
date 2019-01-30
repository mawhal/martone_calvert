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
data <- read.csv( "Output from R/Martone_Hakai_data.csv", stringsAsFactors = FALSE )



# Unique species from the Data
sort( unique( data$Taxon ) )
# write this list of unique names to file
write.csv( data.frame(taxon=sort(unique( data$Taxon ))), "Output from R/Martone_Hakai_uniqueTaxa.csv", row.names=F )

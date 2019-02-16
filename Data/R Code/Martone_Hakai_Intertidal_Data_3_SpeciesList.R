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
data <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data.csv", stringsAsFactors = FALSE )

#metadata
meta<-read.csv("Data/R Code/Output from R/Martone_Hakai_metadata.csv")

#functional group data
functional<-read.csv("Data/taxa/Algae_functional_groups.csv")

# Unique species from the Data
sort( unique( data$Taxon ) )
# write this list of unique names to file
write.csv( data.frame(taxon=sort(unique( data$Taxon ))), "Output from R/Martone_Hakai_uniqueTaxa.csv", row.names=F )

#####Lumping species that are indistinguishable########
#Load lumping data
lump<- read.csv("Data/taxa/CorrectedTaxonList_lumped.csv")

#Create new dataframe with lumped data
lumped.data <- left_join( data, lump, by=c("Taxon"="taxon_corrected") )
algae.lumped <- lumped.data[complete.cases(lumped.data$non.alga.flag=="Algae"),]
algae.lumped2<-left_join(algae.lumped, functional, by=c("Taxon"="Species"))


#Remove non-algae observations from metadata
meta.algae <-meta[complete.cases(lumped.data$non.alga.flag=="Algae"),]
algae.lumped


#Write files
write.csv(algae.lumped2, "Data/R Code/Output from R/Martone_Hakai_data_algae.csv")
write.csv(meta.algae, "Data/R Code/Output from R/Martone_Hakai_metadata_algae.csv")


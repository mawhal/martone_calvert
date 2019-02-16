# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Sam Starko (adapted in part from Matt Whalen "mobr" script)
# updated 16 Feb 2019

# This script is SS's attempt at running community analyses 
# This script uses data on ALGAE ONLY from previous scripts (no animals, no substrate types etc)

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)
library(mobr)
library(vegan)
# data(inv_comm) # Community matrix
# data(inv_plot_attr) # Plot attributes data.frame
# inv_mob_in <- make_mob_in(inv_comm, inv_plot_attr)
# inv_mob_in

## read data files
# all data
algae <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_algae.csv")
# all metadata
algae_metadata <- read.csv("Data/R Code/Output from R/Martone_Hakai_metadata_algae.csv" )

## Data cleaning for Analysis
# remove 2011 data
algae_metadata <- algae_metadata[ algae_metadata$Year != "2011", ]
# remove Meay Channel
algae_metadata <- algae_metadata[ algae_metadata$Site != "Meay Channel", ]

# add up all taxa that appear more than once in a single quadrat (e.g. barnacles or Hildenbrandia) -- go back to datasheets on some of these?
algae[ duplicated(algae), ] # generalize this to look at particular columns

# spread out all of the community data so sample (site,height,year,quadrat) are rows and species are columns
# add together taxa that are not unique to each quadrat
ad.simple <- algae %>%
  group_by( UID, taxon_lumped ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
ad.comm <- ad.simple %>%
  spread( taxon_lumped, Abundance, fill=0 )


am$UID[ !(am$UID %in% ad.comm$UID) ]
ad[ ad$UID %in% am$UID[ !(am$UID %in% ad.comm$UID) ], ]
# restrict to rows selected in metadata
ad.comm <- ad.comm[ ad.comm$UID %in% am$UID, ] 
ad.comm <- as.matrix(ad.comm[,-1])


##Check that the above code resulted in metadata and community files with the same number of rows
length(ad.comm[,1])==length(am[,1])

##Separate community matrices by site
WB.comm <- ad.comm[am$Site=="West Beach",]
WB.meta <- am[am$Site=="West Beach",]
NB.comm <- ad.comm[am$Site=="North Beach",]
NB.meta <- am[am$Site=="North Beach",]
FB.comm <- ad.comm[am$Site=="Fifth Beach",]
FB.meta <- am[am$Site=="Fifth Beach",]

##Diversity indices
##Shannon-Weiner
WB.Shannon <- diversity(WB.comm)
NB.Shannon <- diversity(NB.comm)
FB.Shannon <- diversity(FB.comm)

#Simpson
WB.Simp <- diversity(WB.comm, index="simpson")
NB.Simp <- diversity(NB.comm, index="simpson")
FB.Simp <- diversity(FB.comm, index="simpson")

#Richness
WB.richness <- data.frame()
for (i in 1:length(WB.comm[,1])) {
  sumsum<-sum(WB.comm[i,]>0) 
  WB.richness[i,1]<-sumsum
}

NB.richness <- data.frame()
for (i in 1:length(NB.comm[,1])) {
  sumsum<-sum(NB.comm[i,]>0) 
  NB.richness[i,1]<-sumsum
}


FB.richness <- data.frame()
for (i in 1:length(FB.comm[,1])) {
  sumsum<-sum(FB.comm[i,]>0) 
  FB.richness[i,1]<-sumsum
}

#Combine metadata and new diversity indices
WB.data<-as.data.frame(cbind(WB.meta, WB.Shannon, WB.Simp, WB.richness=as.numeric(unlist(WB.richness))))
NB.data<-cbind(NB.meta, NB.Shannon, NB.Simp, NB.richness=as.numeric(unlist(NB.richness)))
FB.data<-cbind(FB.meta, FB.Shannon, FB.Simp, FB.richness=as.numeric(unlist(FB.richness)))


##Plot species richness versus time
par(mfrow=c(3,1), mar=c(4,4,1,1))
boxplot(WB.richness~Year,data=WB.data[WB.data$Zone=="HIGH",], las=1)
boxplot(WB.richness~Year,data=WB.data[WB.data$Zone=="MID",], las=1)
boxplot(WB.richness~Year,data=WB.data[WB.data$Zone=="LOW",], las=1)

boxplot(NB.richness~Year,data=NB.data[NB.data$Zone=="HIGH",], las=1)
boxplot(NB.richness~Year,data=NB.data[NB.data$Zone=="MID",], las=1)
boxplot(NB.richness~Year,data=NB.data[NB.data$Zone=="LOW",], las=1)

boxplot(FB.richness~Year,data=FB.data[FB.data$Zone=="HIGH",], las=1)
boxplot(FB.richness~Year,data=FB.data[FB.data$Zone=="MID",], las=1)
boxplot(FB.richness~Year,data=FB.data[FB.data$Zone=="LOW",], las=1)


##Plot Shannon versus time
boxplot(WB.Shannon~Year,data=WB.data[WB.data$Zone=="HIGH",], las=1)
boxplot(WB.Shannon~Year,data=WB.data[WB.data$Zone=="MID",], las=1)
boxplot(WB.Shannon~Year,data=WB.data[WB.data$Zone=="LOW",], las=1)

boxplot(NB.Shannon~Year,data=NB.data[NB.data$Zone=="HIGH",], las=1)
boxplot(NB.Shannon~Year,data=NB.data[NB.data$Zone=="MID",], las=1)
boxplot(NB.Shannon~Year,data=NB.data[NB.data$Zone=="LOW",], las=1)

boxplot(FB.Shannon~Year,data=FB.data[FB.data$Zone=="HIGH",], las=1)
boxplot(FB.Shannon~Year,data=FB.data[FB.data$Zone=="MID",], las=1)
boxplot(FB.Shannon~Year,data=FB.data[FB.data$Zone=="LOW",], las=1)

##Plot Simpson versus time
boxplot(WB.Simp~Year,data=WB.data[WB.data$Zone=="HIGH",], las=1)
boxplot(WB.Simp~Year,data=WB.data[WB.data$Zone=="MID",], las=1)
boxplot(WB.Simp~Year,data=WB.data[WB.data$Zone=="LOW",], las=1)

boxplot(NB.Simp~Year,data=NB.data[NB.data$Zone=="HIGH",], las=1)
boxplot(NB.Simp~Year,data=NB.data[NB.data$Zone=="MID",], las=1)
boxplot(NB.Simp~Year,data=NB.data[NB.data$Zone=="LOW",], las=1)

boxplot(FB.Simp~Year,data=FB.data[FB.data$Zone=="HIGH",], las=1)
boxplot(FB.Simp~Year,data=FB.data[FB.data$Zone=="MID",], las=1)
boxplot(FB.Simp~Year,data=FB.data[FB.data$Zone=="LOW",], las=1)

##Rarefaction curves 





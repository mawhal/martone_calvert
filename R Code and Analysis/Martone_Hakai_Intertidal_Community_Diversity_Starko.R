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
library(bayou)

#run summarize functions below before loading Rmisc. 
#calling Rmisc seems to block summarize from doing its job.
library(Rmisc)




## read data files
# all data
algae <- read.csv( "Data/R Code/Output from R/Martone_Hakai_data_algae.csv", stringsAsFactors = FALSE)
# all metadata
algae_metadata <- read.csv("Data/R Code/Output from R/Martone_Hakai_metadata.csv" , stringsAsFactors = FALSE)
###Meta data is currently missing a quadrat from 2016 (Fifth Beach)

#lumped data
lump<-read.csv("Data/R Code/Output from R/Martone_Hakai_data_lumped.csv", stringsAsFactors = FALSE)


## Data cleaning for Analysis
# remove 2011 data
algae_metadata <- algae_metadata[ algae_metadata$Year != "2011", ]
# remove Meay Channel
am <- algae_metadata[ algae_metadata$Site != "Meay Channel", ]


# add up all taxa that appear more than once in a single quadrat (e.g. barnacles or Hildenbrandia) -- go back to datasheets on some of these?
algae[ duplicated(algae), ] # generalize this to look at particular columns


#
#lump.data <- data.frame(cbind(UID=lump$UID, taxon_lumped=lump$taxon_lumped, Abundance=as.numeric(lump$Abundance)))
#lump.data$Abundance<-as.numeric(lump.data$Abundance)

# spread out all of the community data so sample (site,height,year,quadrat) are rows and species are columns
# add together taxa that are not unique to each quadrat
lump.simple <- lump %>%
  group_by( UID, taxon_lumped2 ) %>%
  summarize(Abundance=sum(Abundance,na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
lump.comm <- lump.simple %>%
  spread(taxon_lumped2, Abundance, fill=0 )

#Add in row that has no seaweeds (West Beach High 2016 Q4)
lump.comm<-ungroup(lump.comm)
#lump.comm<-add_row(lump.comm, UID="West Beach HIGH 2016 4")
lump.comm[is.na(lump.comm)] <- 0
lump.comm<-group_by(lump.comm, UID)


am$UID[ !(am$UID %in% lump.comm$UID) ]
lump.comm[ lump.comm$UID %in% am$UID[ !(am$UID %in% lump.comm$UID) ], ]
# restrict to rows selected in metadata
lump.comm <- lump.comm[ lump.comm$UID %in% am$UID, ] 


##Sort metadata and community matrix to be the same order
lump.comm<-lump.comm[order(match(lump.comm$UID, am$UID)),]

##Remove UID column from ad.comm
lump.comm<-as.data.frame(lump.comm[,-1])


##do the same for dataset WITH ONLY algae
#
#ad <- data.frame(cbind(UID=algae$UID, taxon_lumped=algae$taxon_lumped, Abundance=as.numeric(algae$Abundance)))
#ad$Abundance<-as.numeric(ad$Abundance)

# spread out all of the community data so sample (site,height,year,quadrat) are rows and species are columns
# add together taxa that are not unique to each quadrat
ad.simple <- algae %>%
  group_by( UID, taxon_lumped ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T))

# spread Taxon column out into many columns filled with abundance/cover data
ad.comm <- ad.simple %>%
  spread(taxon_lumped, Abundance, fill=0 )

#Add in row that has no seaweeds (West Beach High 2016 Q4)
ad.comm<-ungroup(ad.comm)
ad.comm<-add_row(ad.comm, UID="West Beach HIGH 2016 4")
ad.comm[is.na(ad.comm)] <- 0
ad.comm<-group_by(ad.comm, UID)


am$UID[ !(am$UID %in% ad.comm$UID) ]
ad[ ad$UID %in% am$UID[ !(am$UID %in% ad.comm$UID) ], ]
# restrict to rows selected in metadata
ad.comm <- ad.comm[ ad.comm$UID %in% am$UID, ] 



##Check that the above code resulted in metadata and community 
#files with the same number of rows (both should be 630 or 720 depending on setting above)
str(ad.comm)
str(am)

##Sort metadata and community matrix to be the same order
ad.comm<-ad.comm[order(match(ad.comm$UID, am$UID)),]

##Remove UID column from ad.comm
ad.comm<-as.data.frame(ad.comm[,-1])

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

##West Beach High Zone
WB.high2012 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2012,]
WB.high2013 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2013,]
WB.high2014 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2014,]
WB.high2015 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2015,]
WB.high2016 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2016,]
WB.high2017 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2017,]
WB.high2018 <- WB.comm[WB.meta$Zone=="HIGH"&WB.meta$Year==2018,]


z.high2012 <- specaccum(WB.high2012, method="random", permutations = 1000)
z.high2013 <- specaccum(WB.high2013, method="random", permutations = 1000)
z.high2014 <- specaccum(WB.high2014, method="random", permutations = 1000)
z.high2015 <- specaccum(WB.high2015, method="random", permutations = 1000)
z.high2016 <- specaccum(WB.high2016, method="random", permutations = 1000)
z.high2017 <- specaccum(WB.high2017, method="random", permutations = 1000)
z.high2018 <- specaccum(WB.high2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.high2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,20))
title("High Zone at West Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.high2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.high2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.high2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##West Beach Mid Zone
WB.mid2012 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2012,]
WB.mid2013 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2013,]
WB.mid2014 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2014,]
WB.mid2015 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2015,]
WB.mid2016 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2016,]
WB.mid2017 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2017,]
WB.mid2018 <- WB.comm[WB.meta$Zone=="MID"&WB.meta$Year==2018,]


z.mid2012 <- specaccum(WB.mid2012, method="random", permutations = 1000)
z.mid2013 <- specaccum(WB.mid2013, method="random", permutations = 1000)
z.mid2014 <- specaccum(WB.mid2014, method="random", permutations = 1000)
z.mid2015 <- specaccum(WB.mid2015, method="random", permutations = 1000)
z.mid2016 <- specaccum(WB.mid2016, method="random", permutations = 1000)
z.mid2017 <- specaccum(WB.mid2017, method="random", permutations = 1000)
z.mid2018 <- specaccum(WB.mid2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.mid2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,40))
title("mid Zone at West Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.mid2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.mid2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.mid2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##West Beach Low Zone
WB.low2012 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2012,]
WB.low2013 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2013,]
WB.low2014 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2014,]
WB.low2015 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2015,]
WB.low2016 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2016,]
WB.low2017 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2017,]
WB.low2018 <- WB.comm[WB.meta$Zone=="LOW"&WB.meta$Year==2018,]


z.low2012 <- specaccum(WB.low2012, method="random", permutations = 1000)
z.low2013 <- specaccum(WB.low2013, method="random", permutations = 1000)
z.low2014 <- specaccum(WB.low2014, method="random", permutations = 1000)
z.low2015 <- specaccum(WB.low2015, method="random", permutations = 1000)
z.low2016 <- specaccum(WB.low2016, method="random", permutations = 1000)
z.low2017 <- specaccum(WB.low2017, method="random", permutations = 1000)
z.low2018 <- specaccum(WB.low2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.low2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,60))
title("low Zone at West Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.low2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.low2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.low2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))



##North Beach High Zone
NB.high2012 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2012,]
NB.high2013 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2013,]
NB.high2014 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2014,]
NB.high2015 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2015,]
NB.high2016 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2016,]
NB.high2017 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2017,]
NB.high2018 <- NB.comm[NB.meta$Zone=="HIGH"&NB.meta$Year==2018,]


z.high2012 <- specaccum(NB.high2012, method="random", permutations = 1000)
z.high2013 <- specaccum(NB.high2013, method="random", permutations = 1000)
z.high2014 <- specaccum(NB.high2014, method="random", permutations = 1000)
z.high2015 <- specaccum(NB.high2015, method="random", permutations = 1000)
z.high2016 <- specaccum(NB.high2016, method="random", permutations = 1000)
z.high2017 <- specaccum(NB.high2017, method="random", permutations = 1000)
z.high2018 <- specaccum(NB.high2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.high2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,40))
title("High Zone at North Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.high2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.high2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.high2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##North Beach Mid Zone
NB.mid2012 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2012,]
NB.mid2013 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2013,]
NB.mid2014 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2014,]
NB.mid2015 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2015,]
NB.mid2016 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2016,]
NB.mid2017 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2017,]
NB.mid2018 <- NB.comm[NB.meta$Zone=="MID"&NB.meta$Year==2018,]


z.mid2012 <- specaccum(NB.mid2012, method="random", permutations = 1000)
z.mid2013 <- specaccum(NB.mid2013, method="random", permutations = 1000)
z.mid2014 <- specaccum(NB.mid2014, method="random", permutations = 1000)
z.mid2015 <- specaccum(NB.mid2015, method="random", permutations = 1000)
z.mid2016 <- specaccum(NB.mid2016, method="random", permutations = 1000)
z.mid2017 <- specaccum(NB.mid2017, method="random", permutations = 1000)
z.mid2018 <- specaccum(NB.mid2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.mid2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,50))
title("mid Zone at North Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.mid2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.mid2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.mid2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##North Beach Low Zone
NB.low2011 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2011,]
NB.low2012 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2012,]
NB.low2013 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2013,]
NB.low2014 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2014,]
NB.low2015 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2015,]
NB.low2016 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2016,]
NB.low2017 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2017,]
NB.low2018 <- NB.comm[NB.meta$Zone=="LOW"&NB.meta$Year==2018,]

z.low2011 <- specaccum(NB.low2011, method="random", permutations = 1000)
z.low2012 <- specaccum(NB.low2012, method="random", permutations = 1000)
z.low2013 <- specaccum(NB.low2013, method="random", permutations = 1000)
z.low2014 <- specaccum(NB.low2014, method="random", permutations = 1000)
z.low2015 <- specaccum(NB.low2015, method="random", permutations = 1000)
z.low2016 <- specaccum(NB.low2016, method="random", permutations = 1000)
z.low2017 <- specaccum(NB.low2017, method="random", permutations = 1000)
z.low2018 <- specaccum(NB.low2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.low2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,50))
title("low Zone at North Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.low2011,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.low2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.low2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.low2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))


##Fifth Beach High Zone
FB.high2011 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2011,]
FB.high2012 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2012,]
FB.high2013 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2013,]
FB.high2014 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2014,]
FB.high2015 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2015,]
FB.high2016 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2016,]
FB.high2017 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2017,]
FB.high2018 <- FB.comm[FB.meta$Zone=="HIGH"&FB.meta$Year==2018,]

z.high2011 <- specaccum(FB.high2011, method="random", permutations = 1000)
z.high2012 <- specaccum(FB.high2012, method="random", permutations = 1000)
z.high2013 <- specaccum(FB.high2013, method="random", permutations = 1000)
z.high2014 <- specaccum(FB.high2014, method="random", permutations = 1000)
z.high2015 <- specaccum(FB.high2015, method="random", permutations = 1000)
z.high2016 <- specaccum(FB.high2016, method="random", permutations = 1000)
z.high2017 <- specaccum(FB.high2017, method="random", permutations = 1000)
z.high2018 <- specaccum(FB.high2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.high2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,40))
title("High Zone at Fifth Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.high2011,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.high2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.high2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.high2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.high2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##Fifth Beach Mid Zone
FB.mid2011 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2011,]
FB.mid2012 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2012,]
FB.mid2013 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2013,]
FB.mid2014 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2014,]
FB.mid2015 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2015,]
FB.mid2016 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2016,]
FB.mid2017 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2017,]
FB.mid2018 <- FB.comm[FB.meta$Zone=="MID"&FB.meta$Year==2018,]

z.mid2011 <- specaccum(FB.mid2011, method="random", permutations = 1000)
z.mid2012 <- specaccum(FB.mid2012, method="random", permutations = 1000)
z.mid2013 <- specaccum(FB.mid2013, method="random", permutations = 1000)
z.mid2014 <- specaccum(FB.mid2014, method="random", permutations = 1000)
z.mid2015 <- specaccum(FB.mid2015, method="random", permutations = 1000)
z.mid2016 <- specaccum(FB.mid2016, method="random", permutations = 1000)
z.mid2017 <- specaccum(FB.mid2017, method="random", permutations = 1000)
z.mid2018 <- specaccum(FB.mid2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.mid2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,50))
title("mid Zone at Fifth Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.mid2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.mid2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.mid2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.mid2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

##Fifth Beach Low Zone
FB.low2012 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2012,]
FB.low2013 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2013,]
FB.low2014 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2014,]
FB.low2015 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2015,]
FB.low2016 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2016,]
FB.low2017 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2017,]
FB.low2018 <- FB.comm[FB.meta$Zone=="LOW"&FB.meta$Year==2018,]


z.low2012 <- specaccum(FB.low2012, method="random", permutations = 1000)
z.low2013 <- specaccum(FB.low2013, method="random", permutations = 1000)
z.low2014 <- specaccum(FB.low2014, method="random", permutations = 1000)
z.low2015 <- specaccum(FB.low2015, method="random", permutations = 1000)
z.low2016 <- specaccum(FB.low2016, method="random", permutations = 1000)
z.low2017 <- specaccum(FB.low2017, method="random", permutations = 1000)
z.low2018 <- specaccum(FB.low2018, method="random", permutations = 1000)

par(mfrow=c(1,1), mar=c(5,5,3,3))
plot(z.low2012,ci=1,ci.type="polygon", col="black", lwd=2, ylab="", cex=5,cex.axis=1.1,cex.lab=1.3,las=1, cex.main=1.5,ci.lty=0,ci.col=makeTransparent("blue", 150), ylim=c(0,50))
title("low Zone at Fifth Beach",cex.main=1.75,adj=0,line=1,font.main=7)
title(xlab="Number of Sites",ylab="Richness", line=2, cex.axis=1,cex.lab=1)
plot(z.low2013,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("blue", 200))
plot(z.low2014,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 150))
plot(z.low2015,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2016,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2017,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))
plot(z.low2018,ci=1,add=TRUE, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col=makeTransparent("orange", 200))

####SHORE HEIGHT OF DIFFERENT TRANSECTS
am %>%
  group_by(Zone, Site) %>%
  summarize(mean=mean(Shore_height_cm, na.rm=TRUE), sd=sd(Shore_height_cm, na.rm=TRUE), upper=range(Shore_height_cm, na.rm=TRUE)[1], lower=range(Shore_height_cm)[-1])



###############################################################
#############TOTAL KELP COVER
###############################################################
#############TOTAL KELP COVER
library(Rmisc)

kelps <- c("Alaria marginata", "Costaria costata", "Laminaria setchellii", "Laminaria yezoensis", "Macrocystis pyrifera", "Saccharina sessilis", "Saccharina latissima","Saccharina nigripes", "Egregia menziesii")
kelp.comm<-ad.comm[,kelps]
fucus.comm<-ad.comm[, "Fucus distichus"]
bare.comm<-lump.comm[,"Bare rock"]
pyropia.comm<-lump.comm[,"Pyropia"]
masto.comm<-lump.comm[,"Mastocarpus"]
Ulvales<-c("Ulva", "Blidingia")
ulv.comm<-lump.comm[,Ulvales]
barn.comm<-lump.comm[,"Barnacles"]
head(kelp.comm)
total.kelp <- vector()
for (i in 1:length(kelp.comm[,1])) {
  sumsum<-sum(kelp.comm[i,]) 
  total.kelp[i]<-sumsum
}

total.cover <- vector()
for (i in 1:length(ad.comm[,1])) {
  sumsum<-sum(ad.comm[i,]) 
  total.cover[i]<-sumsum
}

total.ulvales <- vector()
for (i in 1:length(lump.comm[,1])) {
  sumsum<-sum(lump.comm[i,"Ulva"], lump.comm[i,"Blidingia"]) 
  total.ulvales[i]<-sumsum
}

head(total.kelp)
total.kelp<-cbind(total.kelp, am, kelp.comm)
total.cover<-cbind(total.cover, total.kelp)
fucus.cover<-cbind(total.cover, fucus.comm)
bare.rock.cover<-cbind(bare.comm, total.kelp)
pyropia.cover<-cbind(pyropia.comm, total.kelp)
masto.cover<-cbind(masto.comm, total.kelp)
ulvales.cover<-cbind(total.ulvales, total.kelp)
barnacle.cover<-cbind(barn.comm, total.kelp)
kelp.summary<-summarySE(total.kelp, measurevar = "total.kelp", groupvars=c("Year", "Zone", "Site"))
alaria.summary<-summarySE(total.kelp, measurevar = "Alaria marginata", groupvars=c("Year", "Zone", "Site"))
egregia.summary<-summarySE(total.kelp, measurevar = "Egregia menziesii", groupvars=c("Year", "Zone", "Site"))
macrocystis.summary<-summarySE(total.kelp, measurevar = "Macrocystis pyrifera", groupvars=c("Year", "Zone", "Site"))
saccharina.summary<-summarySE(total.kelp, measurevar = "Saccharina sessilis", groupvars=c("Year", "Zone", "Site"))
costaria.summary<-summarySE(total.kelp, measurevar = "Costaria costata", groupvars=c("Year", "Zone", "Site"))
laminaria.yez.summary<-summarySE(total.kelp, measurevar = "Laminaria yezoensis", groupvars=c("Year", "Zone", "Site"))
saccharina.nig.summary<-summarySE(total.kelp, measurevar = "Saccharina nigripes", groupvars=c("Year", "Zone", "Site"))
total.summary<-summarySE(total.cover, measurevar = "total.cover", groupvars=c("Year", "Zone", "Site"))
fucus.summary<-summarySE(fucus.cover, measurevar = "fucus.comm", groupvars=c("Year", "Zone", "Site"))
rock.summary<-summarySE(bare.rock.cover, measurevar="bare.comm", groupvars=c("Year", "Zone", "Site"))
pyropia.summary<-summarySE(pyropia.cover, measurevar="pyropia.comm", groupvars=c("Year", "Zone", "Site"))
masto.summary<-summarySE(masto.cover, measurevar="masto.comm", groupvars=c("Year", "Zone", "Site"))
ulvales.summary<-summarySE(ulvales.cover, measurevar="total.ulvales", groupvars=c("Year", "Zone", "Site"))
barn.summary<-summarySE(barnacle.cover, measurevar="barn.comm", groupvars=c("Year", "Zone", "Site"))

par(mfrow=c(2,1), mar=c(5,5,2,2))
#plot(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="West Beach"&kelp.summary$Zone=="HIGH",], pch=19, ylim=c(0,50), las=1)
plot(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="West Beach"&kelp.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="West Beach"&kelp.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="North Beach"&kelp.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="North Beach"&kelp.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="Fifth Beach"&kelp.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="Fifth Beach"&kelp.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="West Beach"&kelp.summary$Zone=="LOW",], pch=19, ylim=c(0,100), las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="West Beach"&kelp.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="North Beach"&kelp.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="North Beach"&kelp.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="Fifth Beach"&kelp.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(total.kelp~Year,data=kelp.summary[kelp.summary$Site=="Fifth Beach"&kelp.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Total cover
par(mfrow=c(3,1), mar=c(4.1,5,2,1))

plot(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="HIGH",], pch=19, ylim=c(0,200), las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)



plot(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="MID",], pch=19, ylim=c(0,200), las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="LOW",], pch=19, ylim=c(0,200), las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="West Beach"&total.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="North Beach"&total.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(total.cover~Year,data=total.summary[total.summary$Site=="Fifth Beach"&total.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)




##Alaria
plot(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="West Beach"&alaria.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="West Beach"&alaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="North Beach"&alaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="North Beach"&alaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="Fifth Beach"&alaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="Fifth Beach"&alaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="West Beach"&alaria.summary$Zone=="LOW",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="West Beach"&alaria.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="North Beach"&alaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="North Beach"&alaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="Fifth Beach"&alaria.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Alaria marginata`~Year,data=alaria.summary[alaria.summary$Site=="Fifth Beach"&alaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Egregia
plot(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="West Beach"&egregia.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="West Beach"&egregia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="North Beach"&egregia.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="North Beach"&egregia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="Fifth Beach"&egregia.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="Fifth Beach"&egregia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="West Beach"&egregia.summary$Zone=="LOW",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="West Beach"&egregia.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="North Beach"&egregia.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="North Beach"&egregia.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="Fifth Beach"&egregia.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Egregia menziesii`~Year,data=egregia.summary[egregia.summary$Site=="Fifth Beach"&egregia.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Saccharina sessilis (common)
plot(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="HIGH",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="LOW",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="West Beach"&saccharina.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="North Beach"&saccharina.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina sessilis`~Year,data=saccharina.summary[saccharina.summary$Site=="Fifth Beach"&saccharina.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)


##Costaria costata (very rare)
plot(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="West Beach"&costaria.summary$Zone=="MID",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="West Beach"&costaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="North Beach"&costaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="North Beach"&costaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="Fifth Beach"&costaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="Fifth Beach"&costaria.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="West Beach"&costaria.summary$Zone=="LOW",], pch=19, ylim=c(0,200), las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="West Beach"&costaria.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="North Beach"&costaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="North Beach"&costaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="Fifth Beach"&costaria.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Costaria costata`~Year,data=costaria.summary[costaria.summary$Site=="Fifth Beach"&costaria.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Laminaria yezoensis (very rare - WB)
plot(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="West Beach"&laminaria.yez.summary$Zone=="MID",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="West Beach"&laminaria.yez.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="North Beach"&laminaria.yez.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="North Beach"&laminaria.yez.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="Fifth Beach"&laminaria.yez.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="Fifth Beach"&laminaria.yez.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="West Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="West Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="North Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="North Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="Fifth Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Laminaria yezoensis`~Year,data=laminaria.yez.summary[laminaria.yez.summary$Site=="Fifth Beach"&laminaria.yez.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Saccharina nigripes (rare - NB and WB)
plot(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="West Beach"&saccharina.nig.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="West Beach"&saccharina.nig.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="North Beach"&saccharina.nig.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="North Beach"&saccharina.nig.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="Fifth Beach"&saccharina.nig.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="Fifth Beach"&saccharina.nig.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="West Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19, ylim=c(0,30), las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="West Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="North Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="North Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="Fifth Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Saccharina nigripes`~Year,data=saccharina.nig.summary[saccharina.nig.summary$Site=="Fifth Beach"&saccharina.nig.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Macrocystis pyrifera (very rare - only Maey Channel?)
plot(`Macrocystis pyrifera`~Year,data=fucus.summary[macrocystis.summary$Site=="West Beach"&macrocystis.summary$Zone=="MID",], pch=19, ylim=c(0,50), las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="West Beach"&macrocystis.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="North Beach"&macrocystis.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="North Beach"&macrocystis.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="Fifth Beach"&macrocystis.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="Fifth Beach"&macrocystis.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="West Beach"&macrocystis.summary$Zone=="LOW",], pch=19, ylim=c(0,30), las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="West Beach"&macrocystis.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="North Beach"&macrocystis.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="North Beach"&macrocystis.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="Fifth Beach"&macrocystis.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`Macrocystis pyrifera`~Year,data=macrocystis.summary[macrocystis.summary$Site=="Fifth Beach"&macrocystis.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)


##Fucus
plot(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="HIGH",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="MID",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="LOW",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="West Beach"&fucus.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="North Beach"&fucus.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`fucus.comm`~Year,data=fucus.summary[fucus.summary$Site=="Fifth Beach"&fucus.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Bare Rock
plot(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="HIGH",], pch=19, ylim=c(0,100), las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="MID",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="LOW",], pch=19, ylim=c(0,30), las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="West Beach"&rock.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="North Beach"&rock.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`bare.comm`~Year,data=rock.summary[rock.summary$Site=="Fifth Beach"&rock.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Pyropia
plot(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="HIGH",], pch=19, ylim=c(0,30), las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="MID",], pch=19, ylim=c(0,80), las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="LOW",], pch=19, ylim=c(0,30), las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="West Beach"&pyropia.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="North Beach"&pyropia.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`pyropia.comm`~Year,data=pyropia.summary[pyropia.summary$Site=="Fifth Beach"&pyropia.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)


##Mastocarpus
plot(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="HIGH",], pch=19, ylim=c(0,25), las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="MID",], pch=19, ylim=c(0,25), las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="LOW",], pch=19, ylim=c(0,25), las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="West Beach"&masto.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="North Beach"&masto.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`masto.comm`~Year,data=masto.summary[masto.summary$Site=="Fifth Beach"&masto.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Ulvales
plot(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="HIGH",], pch=19, ylim=c(0,10), las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="MID",], pch=19, ylim=c(0,10), las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="LOW",], pch=19, ylim=c(0,10), las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="West Beach"&ulvales.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="North Beach"&ulvales.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`total.ulvales`~Year,data=ulvales.summary[ulvales.summary$Site=="Fifth Beach"&ulvales.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)

##Barnacles
plot(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="HIGH",], pch=19, ylim=c(0,100), las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="HIGH",], pch=19,  las=1, col="blue", lwd=3)


plot(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="MID",], pch=19, ylim=c(0,100), las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="MID",], pch=19,  las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="MID",], pch=19,  las=1, col="blue", lwd=3)

plot(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="LOW",], pch=19, ylim=c(0,100), las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="West Beach"&barn.summary$Zone=="LOW",], pch=19,las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="LOW",], pch=19, las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="North Beach"&barn.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)
points(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="LOW",], pch=19,  las=1, col="blue", cex=2)
lines(`barn.comm`~Year,data=barn.summary[barn.summary$Site=="Fifth Beach"&barn.summary$Zone=="LOW",], pch=19, las=1, col="blue", lwd=3)



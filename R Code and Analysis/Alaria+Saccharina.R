# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 20 December 2018

# This script produces shows patterns of Alaria, Saccharina, and Katharina from on transect
site <- "North Beach"


# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)


## read data files
# all data
ad <- read.csv( "../Data/Excel Files/All years raw/Output from R/Martone_Hakai_data.csv")
# all metadata
am <- read.csv("../Data/Excel Files/All years raw/Output from R/Martone_Hakai_metadata.csv" )


## Deal with trace cover and other oddities
# # replace all commas with periods for Abundance
# ad$Abundance <- as.numeric( gsub( ",", "[.]", ad$Abundance ) )
# change "present" to 0,5
ad$Abundance <- gsub( "present", "0.5", ad$Abundance )
# change 1 of Fucus to "trace"
ad$Abundance <- gsub( "1 on Fucus", "trace", ad$Abundance ) 
# change trace to 0.5% cover
ad$Abundance <- gsub( "t.*", "0.5", ad$Abundance, ignore.case = TRUE ) 

# accept only the first thing if separated by certain characters
# asplit <- strsplit( ad$Abundance, split =c("/|;"))
# ad$Abundance <- unlist( lapply( asplit, function(z) z[1] ))
sort(unique(ad$Abundance))
# things to fix

# spaces
# periods after numbers -- e.g. "0.5."

# taxa to remove...for now
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]


## Identitfy customizations to carry through to figures
# Choose a taxon 
sort( unique( ad$Taxon ))
# use general exp to pull several groups if needed
sacc <- "Saccharina sess"
sort(unique( ad$Taxon[ grep( paste0(sacc,"*"), ad$Taxon ) ]  ))
dsacc   <- ad[ grep( paste0(sacc,"*"), ad$Taxon ), ]
dsacc   <- dsacc[ -grep("Ectocarpus*",dsacc$Taxon), ]

alar <- "Alaria"
sort(unique( ad$Taxon[ grep( paste0(alar,"*"), ad$Taxon ) ]  ))
dalar   <- ad[ grep( paste0(alar,"*"), ad$Taxon ), ]

katy <- "Kat"
sort(unique( ad$Taxon[ grep( paste0(katy,"*"), ad$Taxon ) ]  ))
dkaty   <- ad[ grep( paste0(katy,"*"), ad$Taxon ), ]


# get all instances of a particular taxon
# to include all quadrats us full_join, or use left_join for quads with the taxon
ds <- full_join( dsacc, am, by=c("SiteHeightYear","Quadrat"="Quadrat.No.") )
ds$Abundance[ is.na( ds$Abundance) ] <- 0
ds$group <- sacc
da <- full_join( dalar, am, by=c("SiteHeightYear","Quadrat"="Quadrat.No.") )
da$Abundance[ is.na( da$Abundance) ] <- 0
da$group <- alar
dk <- full_join( dkaty, am, by=c("SiteHeightYear","Quadrat"="Quadrat.No.") )
dk$Abundance[ is.na( dk$Abundance) ] <- 0
dk$group <- katy

# combine all datasets
dask <- do.call( rbind, list(ds,da,dk))

# only look at North Beach, where we have high abundance of all three players
dn <- dask[ dask$Site == site, ]
    
# make abundances numeric
sort(unique(dn$Abundance))
dn$Abundance <- as.numeric( dn$Abundance )
# define leveles for zones
dn$Zone <- factor( dn$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
dn$Site <- factor( dn$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))

# all sites
# time trends in different tidal heights
windows(5,7)
(ggzone <- ggplot( dn, aes(x=Year,y=Abundance)) + facet_grid(group~Zone, scales="free_y") + 
    geom_point(alpha=0.2) + ggtitle( site ) + geom_smooth() )

ggsave( paste0("Figs/AlariaSaccKaty_",site,"_zone.pdf"), ggzone, "pdf" )

# subset of sites where elevation has been measured
windows(10,4)
(ggheight <- ggplot( dn, aes(x=Shore_height_cm,y=Abundance)) + facet_grid(group~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( site ) + geom_smooth(se=FALSE))

ggsave( paste0("Figs/AlariaSaccKaty_",site,"_elevation.pdf"), ggheight, "pdf" )


## correlations

# omit 2011
d11 <- dn[dn$Year>2011,]
# elevation less than 300 (all zeros)
d11 <- d11[d11$Shore_height_cm<300 & !is.na(d11$Shore_height_cm),]


# extract Abundances and reshape
library( reshape2 )
library(psych)
dcor <- d11 %>%
  select( SiteHeightYear, Quadrat, group, Abundance)
dcor <- dcast( dcor, SiteHeightYear + Quadrat ~ group)
dcor <- dcor %>%
  select( Alaria, Kat, Sacc="Saccharina sess")
# get rid of all zero cases for seaweed
dcor <- dcor %>%
  filter( Alaria>0 & Sacc>0 )
pairs.panels(dcor)

ggplot( dcor, aes(x= Alaria, y=Sacc)) + geom_point()

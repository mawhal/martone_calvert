# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# created 07 October 2020

# This script shows patterns of Polysiphonia and Mytilus


# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)


## read data files
# all data
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data.csv" )
# all metadata
am <- read_csv("Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )
# am <- filter(am,Year!=2011)

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
poly <- "Polysiphonia"
sort(unique( ad$Taxon[ grep( paste0(poly,"*"), ad$Taxon ) ]  ))
dpoly   <- ad[ grep( paste0(poly,"*"), ad$Taxon ), ]

myt <- "Mytilus"
sort(unique( ad$Taxon[ grep( paste0(myt,"*"), ad$Taxon ) ]  ))
dmyt   <- ad[ grep( paste0(myt,"*"), ad$Taxon ), ]
dmyt$Taxon <- gsub( " sp.", "", dmyt$Taxon)
dmyt$Taxon <- gsub( " trossulus", "", dmyt$Taxon)


# get all instances of a particular taxon
ds <- left_join( dpoly, am )
ds$Abundance[ is.na( ds$Abundance) ] <- 0
ds$group <- poly
da <- left_join( dmyt, am )
da$Abundance[ is.na( da$Abundance) ] <- 0
da$group <- myt


# combine all datasets
dask <- do.call( rbind, list(ds,da))

# spread and identify transects that had both Mytilus and polyphyllum
dha <- do.call( rbind, list(dpoly,dmyt) ) %>% 
  separate( UID, into = c("beach","x","zone","year","quad"), sep=" ", remove = F ) %>% 
  unite( "transect", beach, x, zone, remove = F ) %>% 
  mutate( Taxon=gsub(" ", "_", Taxon), Abundance=as.numeric(Abundance), year=as.numeric(year)) %>% 
  filter( year != 2011 )
dha_wide <- dha %>% 
  pivot_wider( names_from=Taxon, values_from=Abundance ) 
dha_wide[ is.na(dha_wide) ] <- 0

taxa <- dha_wide[, -c(1:7)]
pairs(taxa)

trans <- dha_wide %>% 
  group_by(transect) %>% 
  summarize( poly=sum(Polysiphonia_hendryi_var._gardneri),myt=sum(Mytilus) ) %>% 
  filter( poly!=0 ) %>% 
  select(transect)
dha_filt <- dha %>% 
  filter( transect %in% trans$transect )
# this contains all transects that had both species and quads with at least one of the species
dha_corr <- dha_filt %>% 
  pivot_wider( names_from=Taxon, values_from=Abundance ) 

# filter out quad with 100% cover of both taxa
dha_corr %>% filter( Mytilus==100 )
dha_corr %>% filter( (Mytilus == 100 & Polysiphonia_hendryi_var._gardneri == 100) ) %>% select(UID)
dha_corr <- dha_corr[ dha_corr$UID != "North Beach LOW 2019 6", ]
dha_corr <- dha_corr %>% 
  mutate( Mytilus = replace_na(Mytilus,replace=0),
          Polysiphonia_hendryi_var._gardneri = replace_na(Polysiphonia_hendryi_var._gardneri,replace=0))

ggplot( dha_corr, aes(x=Mytilus,y=Polysiphonia_hendryi_var._gardneri,col=year) ) + 
  geom_point( alpha=1, size=3 ) +
  # geom_smooth(se=F) +
  facet_wrap(~year) +
  # geom_smooth( aes(group=year), method='lm', se=F ) +
  # stat_ellipse(type="t",level=0.60) +
  theme_bw()
ggsave("R Code and Analysis/Figs/Mytilus+polyphyllum.svg", width=4, height=3)


# correlation, regression
library(lme4)
mm1 <- lmer( Polysiphonia_hendryi_var._gardneri ~ Mytilus + (Mytilus|year/transect), data=dha_corr )
summary(mm1)
    
with( dha_corr, cor.test( Polysiphonia_hendryi_var._gardneri, Mytilus, method="spearman") )


# Cross-correlation
dcf <- dha_corr %>%
  ungroup() %>% 
  select( t=year, Mytilus, Polysiphonia_hendryi_var._gardneri ) 
# is.na(dcf) <- na.omit(dcf)
x <- ccf( dcf$Mytilus, dcf$Polysiphonia_hendryi_var._gardneri, ylab="cross-correlation", main="Mytilus & polyphyllum", lag.max = 60)

#























# make abundances numeric
sort(unique(dn$Abundance))
dn$Abundance <- as.numeric( dn$Abundance )
# define leveles for zones
dn$Zone <- factor( dn$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
dn$Site <- factor( dn$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))
# dn <- filter(dn,!is.na(Year))

# all plots with Mytilus or Polysiphoniaile
dah <- dn %>% filter(Taxon %in% c("Mytilus marginata","Polysiphoniaile")) %>% 
  select( Site, Zone, Year, Meter.point ) %>% distinct
# for every year and each taxon
bits <- expand.grid( Year=2011:2019, Taxon=c("Mytilus marginata","Polysiphoniaile", "Katharina tunicata") )
dnfill <- dn %>% group_by(Year) %>% summarize(Date = mean(lubridate::ymd(Date),na.rm=T))
bits <- full_join(dnfill, bits)
filling <- merge( dah, bits )

# all sites
# time trends in different tidal heights
windows(6,2)
# fl16 <- filter(dn, Taxon %in% c("Mytilus marginata","Polysiphoniaile"), Meter.point==16, Site=="Fifth Beach", Zone=="LOW")
# fl16$Date <- lubridate::ymd(fl16$Date)
toplot <- left_join(filling,dn, by=c("Site","Zone","Meter.point","Year","Taxon"))
toplot$Abundance[is.na(toplot$Abundance)] <- 0
toplot$Abundance[toplot$Taxon =="Katharina tunicata" & toplot$Year==2011] <- NA

ggplot( toplot, aes(x=Date.x,y=Abundance)) + facet_wrap(~Taxon, scales="free_y") + 
    geom_point(alpha=0.2) + geom_smooth() + 
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_classic() +
  theme(axis.title.x=element_blank())
ggsave( "R Code and Analysis/Figs/MytiluspolyKaty.svg", width=6, height=2 )

# subset of sites where elevation has been measured
windows(10,4)
(ggheight <- ggplot( dn, aes(x=Shore_height_cm,y=Abundance)) + facet_grid(group~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( site ) + geom_smooth(se=FALSE))

ggsave( paste0("Figs/MytiluspolyKaty_",site,"_elevation.pdf"), ggheight, "pdf" )


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
  select( Mytilus, Kat, poly="Polysiphonia")
# get rid of all zero cases for seaweed
dcor <- dcor %>%
  filter( Mytilus>0 & poly>0 )
pairs.panels(dcor)

ggplot( dcor, aes(x= Mytilus, y=poly)) + geom_point()

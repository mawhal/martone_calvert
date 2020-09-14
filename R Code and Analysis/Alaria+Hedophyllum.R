# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 20 July 2020

# This script produces shows patterns of Alaria, Hedophyllum sessile, and Katharina from on transect
# site <- "North Beach"


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
hedo <- "Hedophyllum sess"
sort(unique( ad$Taxon[ grep( paste0(hedo,"*"), ad$Taxon ) ]  ))
dhedo   <- ad[ grep( paste0(hedo,"*"), ad$Taxon ), ]

alar <- "Alaria"
sort(unique( ad$Taxon[ grep( paste0(alar,"*"), ad$Taxon ) ]  ))
dalar   <- ad[ grep( paste0(alar,"*"), ad$Taxon ), ]

katy <- "Kat"
sort(unique( ad$Taxon[ grep( paste0(katy,"*"), ad$Taxon ) ]  ))
dkaty   <- ad[ grep( paste0(katy,"*"), ad$Taxon ), ]


# get all instances of a particular taxon
ds <- left_join( dhedo, am )
ds$Abundance[ is.na( ds$Abundance) ] <- 0
ds$group <- hedo
da <- left_join( dalar, am )
da$Abundance[ is.na( da$Abundance) ] <- 0
da$group <- alar
dk <- left_join( dkaty, am )
dk$Abundance[ is.na( dk$Abundance) ] <- 0
dk$group <- katy

# combine all datasets
dask <- do.call( rbind, list(ds,da,dk))

# spread and identify transects that had both Alaria and Hedophyllum
dha <- do.call( rbind, list(dhedo,dalar) ) %>% 
  separate( UID, into = c("beach","x","zone","year","quad"), sep=" ", remove = F ) %>% 
  unite( "transect", beach, x, zone, remove = F ) %>% 
  mutate( Taxon=gsub(" ", "_", Taxon), Abundance=as.numeric(Abundance), year=as.numeric(year)) %>% 
  filter( year != 2011 )
dha_wide <- dha %>% 
  pivot_wider( names_from=Taxon, values_from=Abundance ) 
dha_wide[ is.na(dha_wide) ] <- 0
trans <- dha_wide %>% 
  group_by(transect) %>% 
  summarize( Hedo=sum(Hedophyllum_sessile),Alar=sum(Alaria_marginata) ) %>% 
  filter( Hedo!=0 ) %>% 
  select(transect)
dha_filt <- dha %>% 
  filter( transect %in% trans$transect )
# this contains all transects that had both species and quads with at least one of the species
dha_corr <- dha_filt %>% 
  pivot_wider( names_from=Taxon, values_from=Abundance ) 

# filter out quad with 100% cover of both taxa
dha_corr %>% filter( Alaria_marginata==100 )
dha_corr %>% filter( (Alaria_marginata == 100 & Hedophyllum_sessile == 100) ) %>% select(UID)
dha_corr <- dha_corr[ dha_corr$UID != "North Beach LOW 2019 6", ]
dha_corr <- dha_corr %>% 
  mutate( Alaria_marginata = replace_na(Alaria_marginata,replace=0),
          Hedophyllum_sessile = replace_na(Hedophyllum_sessile,replace=0))

ggplot( dha_corr, aes(x=Alaria_marginata,y=Hedophyllum_sessile,col=year) ) + 
  geom_point( alpha=1, size=3 ) +
  # geom_smooth(se=F) +
  facet_wrap(~year) +
  # geom_smooth( aes(group=year), method='lm', se=F ) +
  # stat_ellipse(type="t",level=0.60) +
  theme_bw()
ggsave("R Code and Analysis/Figs/Alaria+Hedophyllum.svg", width=4, height=3)


# correlation, regression
library(lme4)
mm1 <- lmer( Hedophyllum_sessile ~ Alaria_marginata + (Alaria_marginata|year/transect), data=dha_corr )
summary(mm1)
    
with( dha_corr, cor.test( Hedophyllum_sessile, Alaria_marginata, method="spearman") )


# Cross-correlation
dcf <- dha_corr %>%
  ungroup() %>% 
  select( t=year, Alaria_marginata, Hedophyllum_sessile ) 
# is.na(dcf) <- na.omit(dcf)
x <- ccf( dcf$Alaria_marginata, dcf$Hedophyllum_sessile, ylab="cross-correlation", main="Alaria & Hedophyllum", lag.max = 60)

#























# make abundances numeric
sort(unique(dn$Abundance))
dn$Abundance <- as.numeric( dn$Abundance )
# define leveles for zones
dn$Zone <- factor( dn$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
dn$Site <- factor( dn$Site, levels = c("Fifth Beach", "West Beach", "North Beach", "Meay Channel"))
# dn <- filter(dn,!is.na(Year))

# all plots with Alaria or Hedophyllum sessile
dah <- dn %>% filter(Taxon %in% c("Alaria marginata","Hedophyllum sessile")) %>% 
  select( Site, Zone, Year, Meter.point ) %>% distinct
# for every year and each taxon
bits <- expand.grid( Year=2011:2019, Taxon=c("Alaria marginata","Hedophyllum sessile", "Katharina tunicata") )
dnfill <- dn %>% group_by(Year) %>% summarize(Date = mean(lubridate::ymd(Date),na.rm=T))
bits <- full_join(dnfill, bits)
filling <- merge( dah, bits )

# all sites
# time trends in different tidal heights
windows(6,2)
# fl16 <- filter(dn, Taxon %in% c("Alaria marginata","Hedophyllum sessile"), Meter.point==16, Site=="Fifth Beach", Zone=="LOW")
# fl16$Date <- lubridate::ymd(fl16$Date)
toplot <- left_join(filling,dn, by=c("Site","Zone","Meter.point","Year","Taxon"))
toplot$Abundance[is.na(toplot$Abundance)] <- 0
toplot$Abundance[toplot$Taxon =="Katharina tunicata" & toplot$Year==2011] <- NA

ggplot( toplot, aes(x=Date.x,y=Abundance)) + facet_wrap(~Taxon, scales="free_y") + 
    geom_point(alpha=0.2) + geom_smooth() + 
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_classic() +
  theme(axis.title.x=element_blank())
ggsave( "R Code and Analysis/Figs/AlariahedoKaty.svg", width=6, height=2 )

# subset of sites where elevation has been measured
windows(10,4)
(ggheight <- ggplot( dn, aes(x=Shore_height_cm,y=Abundance)) + facet_grid(group~Year, scales = "free_y") + 
    geom_point(alpha=0.2) +  ggtitle( site ) + geom_smooth(se=FALSE))

ggsave( paste0("Figs/AlariahedoKaty_",site,"_elevation.pdf"), ggheight, "pdf" )


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
  select( Alaria, Kat, Hedo="Hedophyllum sess")
# get rid of all zero cases for seaweed
dcor <- dcor %>%
  filter( Alaria>0 & Hedo>0 )
pairs.panels(dcor)

ggplot( dcor, aes(x= Alaria, y=Hedo)) + geom_point()

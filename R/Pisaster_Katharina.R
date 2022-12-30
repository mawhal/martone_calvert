# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# This script produces figures of the density of a chosen taxa, saving figures as pdf
taxon <- "Pisaster"
taxon2 <- "Katharina"

# set options
options(stringsAsFactors = FALSE)

# load packages
require(svglite)
library(tidyverse)
library(viridis)

## read data files
# all data
ad <- read.csv( "Data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv")
sort(unique( ad$taxon_lumped3[ grep( paste0(taxon,"*"), ad$taxon_lumped3 ) ]  ))
sort(unique( ad$taxon_lumped3[ grep( paste0(taxon2,"*"), ad$taxon_lumped3 ) ]  ))
# all metadata
am <- read.csv("Data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )




# taxa to remove
ad <- ad[ ad$Taxon != "Black spots on Fucus", ]

# customize sites and years
years <- 2012:2019
sites <- c("Fifth Beach", "North Beach","Foggy Cove")
metause <- am %>% 
  filter( Year %in% years ) %>% 
  filter( Site %in% sites )

## Customizations to carry through to figures

# Choose a taxon 
# use general exp to pull several groups if needed
# dtax   <- ad[ grep( paste0(taxon,"*"), ad$taxon_revised ), ]
dtax   <- ad[ grep( paste0(taxon,"*"), ad$taxon_lumped3 ), ]
dtax2   <- ad[ grep( paste0(taxon2,"*"), ad$taxon_lumped3 ), ]


# merge data and metadata
# get all instances of a particular taxon, and code missing as zero
# include all quadrats
d <- left_join( metause, dtax )
d2 <- left_join( metause, dtax2 )
dd <- bind_rows( d, d2 )
d <- dd
# # get rid of bad merging
# d <- d %>% 
#   filter( !is.na(Site), !is.na(Zone) )
# # query d for consistent plots revisited
# drevisited <- d %>% 
#   group_by( Site, Zone, Meter.point ) %>% 
#   summarize( n = length(Quadrat) ) %>% 
#   arrange( -n )
# dfilt <- left_join(d, drevisited) %>% 
#   filter( n > 5 )
# d <- dfilt

# *** to add *** include full names for each unique instance of Sampler


# chose an observer/recorder
sort(unique( d$Sampler ))
sort(unique( d$Recorder ))
allpeeps <- data.frame( Person = sort(unique( c(d$Sampler,d$Recorder))), fullName=""  )
dsampler <- d %>% 
  group_by( Year, Sampler, Recorder) %>% 
  summarize( samples = n() )
dperson <- dsampler %>% 
  pivot_longer( cols = c("Sampler","Recorder"), names_to = "Person", values_to = "sample" )
dperson <- dperson %>% 
  select( person = sample, Year, samples ) %>% 
  group_by( Year, person ) %>% 
  summarize( samples = sum(samples) )


dperson$person[ dperson$person %in% c('SL','SCL?') ] <- 'SCL'
sort(unique(dperson$person))
with(dperson, table(person, Year) )


# Sandra C Lindstrom
sandra     <- c( "Sandra", "SCL", "SL" )
# select the plots we want from the metadata
dperson <- d %>%
  filter( Sampler %in% sandra | Recorder %in% sandra )
d<- dperson





# make abundances numeric
sort(unique(d$Abundance))
d$Abundance <- as.numeric( d$Abundance )
# define leveles for zones
d$Zone <- factor( d$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
d$Site <- factor( d$Site, levels = sites )
# define Year factor
d$Year <- factor( d$Year, ordered= TRUE )
#Remove taxa that are NA
d=d[!is.na(d$Taxon),]

# completely fill all observations (adding zeros to quadrats where available)
dplot <- d %>% 
  select( Year, Site, Zone, Quadrat, Taxon, Abundance)# %>% 
  # complete( Year, Site, Zone, Quadrat, fill = list(Abundance = 0) ) %>% 
  # filter( !(Year %in% 2011:2014) | Site != "Meay Channel" )



# make a numeric year
dplot$yearnum <- as.numeric(as.character(dplot$Year))

# all sites
# time trends in different tidal zones
# quartz(5,5)

(ggzone <- ggplot( dplot, aes(x=yearnum,y=Abundance, col = Taxon)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( aes(group=Taxon), se=TRUE, col='black' ) +
    stat_summary( aes(group=Taxon), fun.data = "mean_cl_boot",  geom = "errorbar", size = 0.5, width = 0.2 ) +
    stat_summary( aes(group=Taxon), fun = "mean", geom="line", size = 0.5 ) +
    geom_point( alpha=0.4 ) + 
    xlab("Year") + ylab("Count") +
    scale_x_continuous(breaks = seq(2010,2022,2) ) )


# widen data for correlation
dcorr <- dplot %>% 
  pivot_wider( names_from = Taxon, values_from = Abundance )
dcorr <- dcorr[,7:8]
dcorr <- dcorr[ complete.cases(dcorr), ]
cor(dcorr)
cor.test( dcorr$`Pisaster ochraceus`, dcorr$`Katharina tunicata`)
# ggsave( paste0("R/Figs/",taxon,"_zone.svg") )

# plot means across the dataset for transect in which the taxon of interest appeared at least once
# pick transects
dplot <- dplot %>% 
  unite( transect, Site, Zone, remove=F )
# summarize counts by quadrat
dquad <- dplot  %>% 
  group_by( Year, Taxon, yearnum, Site, Zone, transect, Quadrat ) %>% 
  summarize( Abundance=sum(Abundance) ) 

(ggzone <- ggplot( dquad, aes(x=yearnum,y=Abundance, col=Taxon)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=TRUE, col='black' ) +
    stat_summary( fun.data = "mean_cl_boot",  geom = "errorbar", colour = "slateblue4", size = 0.5, width = 0.2 ) +
    stat_summary( fun = "mean", geom="line", colour = "slateblue4", size = 0.5 ) +
    geom_point( alpha=0.4,col='slateblue' ) + 
    xlab("Year") + ylab("Count") +
    scale_x_continuous(limits = c(2012,2019), 
                       breaks = seq(2010,2022,2) ) )


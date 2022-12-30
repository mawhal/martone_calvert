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


# # chose an observer/recorder
# sort(unique( am$Sampler ))
# sort(unique( am$Recorder ))
# allpeeps <- data.frame( Person = sort(unique( c(am$Sampler,am$Recorder))), fullName=""  )
# write_csv(allpeeps, "R/output/all_samplers_recorders.csv")
# # Sandra C Lindstrom
# sandra     <- c( "Sandra", "SCL", "SL" )
# # select the plots we want from the metadata
# dperson <- d %>%
#   filter( Sampler %in% sandra | Recorder %in% sandra )
# d<- dperson





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
  select( Year, Site, Zone, Quadrat, Abundance) %>% 
  # complete( Year, Site, Zone, Quadrat, fill = list(Abundance = 0) ) %>% 
  filter( !(Year %in% 2011:2014) | Site != "Meay Channel" )



# make a numeric year
dplot$yearnum <- as.numeric(as.character(dplot$Year))

# all sites
# time trends in different tidal zones
# quartz(5,5)

(ggzone <- ggplot( dplot, aes(x=yearnum,y=Abundance)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=TRUE, col='black' ) +
    stat_summary( fun.data = "mean_cl_boot",  geom = "errorbar", colour = "slateblue4", size = 0.5, width = 0.2 ) +
    stat_summary( fun = "mean", geom="line", colour = "slateblue4", size = 0.5 ) +
    geom_point( alpha=0.4,col='slateblue' ) + ggtitle( taxon ) + 
    xlab("Year") + ylab("Count") +
    scale_x_continuous(breaks = seq(2010,2022,2) ) )

# ggsave( paste0("R/Figs/",taxon,"_zone.svg") )

# plot means across the dataset for transect in which the taxon of interest appeared at least once
# pick transects
dplot <- dplot %>% 
  unite( transect, Site, Zone, remove=F )
# summarize counts by quadrat
dquad <- dplot  %>% 
  group_by( Year, yearnum, Site, Zone, transect, Quadrat ) %>% 
  summarize( Abundance=sum(Abundance) ) 

(ggzone <- ggplot( dquad, aes(x=yearnum,y=Abundance)) + 
    facet_grid(Site~Zone, scales="free_y") + 
    # geom_smooth( se=TRUE, col='black' ) +
    stat_summary( fun.data = "mean_cl_boot",  geom = "errorbar", colour = "slateblue4", size = 0.5, width = 0.2 ) +
    stat_summary( fun = "mean", geom="line", colour = "slateblue4", size = 0.5 ) +
    geom_point( alpha=0.4,col='slateblue' ) + ggtitle( taxon ) + 
    xlab("Year") + ylab("Count") +
    scale_x_continuous(limits = c(2012,2019), 
                       breaks = seq(2010,2022,2) ) )


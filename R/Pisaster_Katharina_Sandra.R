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

##
# chose an observer/recorder
sort(unique( metause$Sampler ))
sort(unique( metause$Recorder ))
allpeeps <- data.frame( Person = sort(unique( c(metause$Sampler,metause$Recorder))), fullName=""  )
# dsampler <- d %>% 
#   group_by( Year, Sampler, Recorder) %>% 
#   summarize( samples = n() )
# dperson <- dsampler %>% 
#   pivot_longer( cols = c("Sampler","Recorder"), names_to = "Person", values_to = "sample" )
# dperson <- dperson %>% 
#   select( person = sample, Year, samples ) %>% 
#   group_by( Year, person ) %>% 
#   summarize( samples = sum(samples) )
sort(unique(allpeeps$Person))
# Sandra C Lindstrom
sandra     <- c( "Sandra", "SCL", "SL", "SCL?" )
# select the plots we want from the metadata
metause <- metause %>% 
  filter( Sampler %in% sandra | Recorder %in% sandra )


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
d <- d %>% select( Site, Zone, Year, Quadrat, Taxon, Abundance )
d$Taxon[ is.na(d$Taxon) ] <- "Pisaster ochraceus"
d$Abundance[ is.na(d$Abundance) ] <- 0
d2 <- left_join( metause, dtax2 )
d2 <- d2 %>% select( Site, Zone, Year, Quadrat, Taxon, Abundance )
d2$Taxon[ is.na(d2$Taxon) ] <- "Katharina tunicata"
d2$Abundance[ is.na(d2$Abundance) ] <- 0
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


# calculate sample size
sandra.samples <- metause %>% 
  group_by( Year, Site, Zone ) %>% 
  summarize( sample.size = n() ) %>% 
  mutate( yearnum = as.numeric(Year) )

# define leveles for zones
dplot$Zone <- factor( dplot$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
sandra.samples$Zone <- factor( sandra.samples$Zone, levels = c("LOW","MID","HIGH"), ordered = T )
# define Site order
dplot$Site <- factor( dplot$Site, levels = sites )

sandra.samples$Taxon <- "Pisaster ochraceus"
sandra.samples <- bind_rows( sandra.samples, data.frame( Year = 2014, Site = "North Beach", Zone = "MID", sample.size = 0, yearnum = 2014, Taxon = "Pisaster ochraceus" ))
sandra.samples$Site <- factor( sandra.samples$Site, levels = sites )
sandra.samples$Zone <- factor( sandra.samples$Zone, levels = c("LOW","MID","HIGH"), ordered = T )

(ggzone <- ggplot( dplot, aes(x=yearnum,y=Abundance, col = Taxon, shape = Taxon)) + 
    facet_grid(Site~Zone) + 
    stat_summary( aes(group=Taxon), fun = "mean", geom="line", size = 0.5 ) +
    geom_point( alpha=0.33 ) + 
    geom_text( data = sandra.samples, aes( x = yearnum, y = 6, label = sample.size),
               col = "black", size = 2.5 ) +
    xlab("Year") + ylab("Count") +
    scale_x_continuous(breaks = seq(2010,2022,2) ) ) +
  theme( legend.position = 'top', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )
ggsave( "R/Figs/pisaster_katharina_sandra.svg", width = 4, height = 4 )


# widen data for correlation
dwide <- dplot %>% 
  pivot_wider( names_from = Taxon, values_from = Abundance )
dcorr <- dwide[,6:7]
dcorr <- dcorr[ complete.cases(dcorr), ]
cor(dcorr)
cor.test( dcorr$`Pisaster ochraceus`, dcorr$`Katharina tunicata`)
cor.test( dcorr$`Pisaster ochraceus`, dcorr$`Katharina tunicata`, method = "spearman")
plot(dcorr)
ggplot( dwide, aes(x = `Pisaster ochraceus`, y = `Katharina tunicata`)) + 
  geom_point(alpha = 0.25)
psych::pairs.panels(dcorr)


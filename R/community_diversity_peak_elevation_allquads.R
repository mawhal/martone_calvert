# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen 
# created 2 Dec 2021

# Using Wang et al. 2019 framework for stability and synchrony
# https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.04290

# 

# 
library(tidyverse)
# library(codyn)
library(str2str)
## var.partition function
# source( "R/var.partition.R" )

## Evenness as defined as Evar in Smith & Wilson 1996 Oikos
Evar <- function( x ){
  S = length( x[x>0] )
  1 - 2/pi*atan( sum((log(x[ x>0 ]) - sum(log(x[ x>0 ]))/S)^2)/S ) 
}


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
am <- am[ am$Year != "2011", ]
# remove Meay Channel
am <- am[ am$Site != "Meay Channel", ]

# read in data from Community_rda.R
d.simple <- read_csv("R/output/data_select_rda_HMSC.csv")

# define quadrats
d.simple <- d.simple %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point) )


# just get algae
d.comm.algae <- d.simple %>%
  filter( funct_2021 != "animal" ) %>% 
  select(-funct_2021) %>% 
  spread( taxon, Abundance, fill=0 )





# merge meta data so we can chop things up and summarize across sites, zones, etc.
# first, remove rows from data that are not in the restricted metadata
am.select  <- am %>% 
  mutate( quadrat = paste( Site, Zone, Meter.point, sep = " ") ) %>% 
  select( quadrat, Site, Zone, Meter.point, Shore_height_cm )
muse <- left_join( d.simple, distinct(am.select) )
splits <- strsplit( as.character(muse$UID), " " )
muse$transect <- unlist(lapply(splits, function(z) paste(z[1:4],collapse = " ")))



# # there is a quadrat without any metadata, remove this from the metadata
# noalgae <- anti_join( muse, d.comm )
# noalgae$UID
# mclean <- muse[ muse$UID != noalgae$UID, ]
mclean <- muse
# define levels for zones
mclean$Zone <- factor( mclean$Zone, levels = c("LOW","MID","HIGH"), ordered = F )
# define Site order
mclean$Site <- factor( mclean$Site, levels = c("Foggy Cove", "Fifth Beach", "North Beach" ))
# define transects
mclean <- mclean %>% 
  unite( "transect", Site, Zone, remove = F, sep = " ")

##Sort metadata and community matrix to be the same order
# d.comm.order <- d.comm[ order(match(d.comm$transect, mtrans$transect)),]
# cbind( d.comm.order$transect, mtrans$transect )

# remove UID column from community data
# comm.all <- as.matrix(d.comm.all[,-c(1:6)])
comm.algae <- as.matrix(d.comm.algae[,-c(1:7)])





melev <-  mclean %>% 
  select( UID, Year, transect, Site, Zone, quadrat, Shore_height_cm ) %>% 
  distinct()
  
dd <- left_join(melev, d.comm.algae)




# remove taxa that do not appear
colSums(dd[,-c(1:9)])
ddmean.transect <- dd %>% 
  separate( quadrat, c('Site', 'blah', 'Zone', 'meter.point'), sep=" ", remove = F ) %>% 
  unite( "transect", Site, blah, Zone, sep = " ") %>% 
  mutate( transect = factor(transect) )









# Now come back to HMSC analysis. Let's consider the peak elevation of a species at the start of the survey (first two years?)
elev_abun_shifts <- read_csv( "R/output/shifts_predicted.csv" )

# compare variability of elevation optima within transects to stability and synchrony
# replace occurrences with initial elevation peaks
# variability as a mean-variance relationship?
initial_peaks_taxon <- elev_abun_shifts %>% 
  select( taxon, elev.init.med ) 



# compare variability using composition data from first two years
cover_quads_all <- ddmean.transect 

cover_quads_all_comm <- select( cover_quads_all, Acrosiphonia:Unknown.crust )
cover_quads_all_pa <- ifelse( cover_quads_all_comm > 0, 1, 0)

# only select columns for which we have results from HMSC
pa_hmsc <- cover_quads_all_pa[, colnames(cover_quads_all_pa) %in% initial_peaks_taxon$taxon ]
# filter out peak elevations missing from surveyed quadrats
# interesting that we are missing 3 taxa from 85 quadrats sampled over four years...but, goes to show ya, eh?
initial_peaks_taxon <- filter( initial_peaks_taxon, taxon %in% colnames(pa_hmsc) )

# get things in same order
initial_peaks_taxon <- initial_peaks_taxon %>% 
  arrange(taxon)

# for loop
all( initial_peaks_taxon$taxon == colnames(pa_hmsc) )
for(i in 1:ncol(pa_hmsc) ){
  pa_hmsc[,i] <-  ifelse( pa_hmsc[,i] == 1, initial_peaks_taxon$elev.init.med[i], NA )
}

# calculate mean and sd of peak elevation for each quadrat
peaks.mean <- apply( pa_hmsc, 1, mean, na.rm = T)
peaks.sd <- apply( pa_hmsc, 1, sd, na.rm = T)
boxplot(peaks.sd/peaks.mean)




cover_quads_all$elev.var <- peaks.sd/peaks.mean
cover_quads_all$localrich <- rowSums(cover_quads_all_pa)

# windows( 5,5 )
cover_quads_pre <- filter( cover_quads_all, Year %in% c(2012,2013))
psych::pairs.panels( select(cover_quads_pre, Shore_height_cm, localrich, elev.var), hist.col = "whitesmoke"                )




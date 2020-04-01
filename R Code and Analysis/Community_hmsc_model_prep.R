# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space


# load libraries
library( tidyverse )
library( Hmsc )
library( corrplot )
library( here )
library( tictoc )

# useful references for this code
citation( 'Hmsc' )


# dataDir = file.path(localDir, "data")
ModelDir = paste0( here::here(), "/R Code and Analysis/models" )
MixingDir = paste0( here::here(), "/R Code and Analysis/mixing")



## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read.csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv", stringsAsFactors = FALSE )
# all metadata
am <- read.csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )





## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
# muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]

# no NA values allowed, so we need to remove these from the dataset
rem <- unique( which(is.na(muse$Shore_height_cm))  )
# get rid of row 2 in all data structures
muse   <- muse[-rem,]


muse <- droplevels(muse)


# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )


# for now, restrict community analysis to sessile organisms (mobile data not very reliable)
d <- dm %>% 
  filter( motile_sessile != "motile" )
  # filter( non.alga.flag =="Algae" )



# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  group_by( UID, Year, Site, Zone, taxon_lumped2 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) %>%
  mutate( taxon_lumped2 = gsub(" ",".",taxon_lumped2) )

# # average cover per transect
# dmean <- d.simple %>%
#   spread( taxon_lumped2, Abundance, fill=0 ) %>%
#   gather( taxon_lumped2, Abundance, -UID, -Year, -Site, -Zone ) %>%
#   group_by( Year, Site, Zone, taxon_lumped2 ) %>%
#   summarise( Abundance=mean(Abundance) )

d.comm.prep <- d.simple  # dmean

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.comm.prep %>%
  spread( taxon_lumped2, Abundance, fill=0 )

# order community data by site and zone
d.comm <- d.comm %>%
  arrange( Site, factor(Zone,levels=c("LOW","MID","HIGH")) )
d.comm$Zone <- factor( d.comm$Zone,levels=c("LOW","MID","HIGH") )



# isolate the community, site, and sample data
comm.all <- d.comm[,-c(1:4)]
# colSums(comm.all)
# summary(colSums(comm.all))
# hist(colSums(comm.all))
# sort(colSums(comm.all))
# sort(colSums(comm.all),decreasing = TRUE)[1:10]
# boxplot( colSums(comm.all))
# boxplot( log(colSums(comm.all)) )
# boxplot( colSums(comm.all)[ colSums(comm.all)> 10] ) 
# hist( colSums(comm.all)[ colSums(comm.all)> 10] )
# sort(colSums(comm.all)[ which( colSums(comm.all) <=10 )  ])
# length(colSums(comm.all)[ which( colSums(comm.all) <=10 )  ])

# reduce the dataset by removing the rarest taxa 
# those that have less than a total percent cover 
# acheived across all of the selected quadrats 
# and across all selected years
which( colSums(comm.all) <=10 )
# filter by occurence
occurrence <- apply(comm.all, 2, function(z) length(z[z>0])) 
comm <- comm.all %>% 
  select( names(occurrence)[ which(occurrence >= length(unique(d.comm$Year))*6) ] )
# filter by abundance
abundance <- colSums(comm)d
which( colSums(comm) <=10 )
# comm <- comm[ ,-which( abundance <=10 )  ]
# comm <- comm.all
commpa <- comm
commpa[commpa>0] <- 1
sort(colSums(commpa))[1:20]
# reorder community matrix based on the most abundant to least (across all quadrats)
Y <- as.matrix( comm[, order(colSums(comm),decreasing = T) ] )

# compare Y and ad
ncol(Y)
ad %>% filter(motile_sessile == "sessile") %>% 
  select(taxon_lumped2,motile_sessile,kelp_fucoid_turf,new_cat_simple3) %>%
  arrange( taxon_lumped2 ) %>% 
  distinct()

ad %>% filter(motile_sessile == "sessile") %>% 
  select(taxon_lumped2) %>%
  arrange( taxon_lumped2 ) %>% 
  distinct() #%>% 
ad.rich <- ad %>% filter(motile_sessile == "sessile") %>% 
  select(taxon_lumped2) %>%
  arrange( taxon_lumped2 ) %>% 
  distinct() %>% 
  summarize(length(unique(taxon_lumped2)))

ncol(Y)/ad.rich


# define the variables to test from metadata and data
# merge community data with metadata
metacomm <- left_join( d.comm, muse )
write_csv( metacomm, "R Code and Analysis/output from r/community.csv")

# read temperature data
pine <- read_csv( "R Code and Analysis/output from r/PineIsland_summary.csv" )

# merge temperature and the rest of the metadata
M <- left_join( metacomm, pine )
# qudrat as another potential term? but should be captured by lat,long
M$transect <- unlist(lapply( strsplit( M$UID, " " ), function(z) paste( z[1],z[2],z[3] ) ))

# # random?
# year
# geolocation
# 
# # fixed?
# year
# temperature
# shore height
# 

## Fixed Effects
XData <- M %>%
  ungroup() %>%
  select( year=Year, shore.height = Shore_height_cm, anom.pine.sum.1, anom.pine.win ) #%>% 
  # mutate( year=factor(year, ordered=TRUE ) )
XData <- as.data.frame( XData )


# set up random effects for hmsc
# make sure to use spatial data
# STUDY DESIGN
studyDesign <- data.frame( year = factor(M$Year),
                           observation = factor(M$UID),
                           transect = factor(M$transect),
                           site = factor(M$Site) )
# transect year
studyDesign$ty <- factor(with(studyDesign, paste( transect, year, sep = "_" )))


# Random effects for unit
# rL_obs <- HmscRandomLevel( unit = studyDesign$observation )
rL <- HmscRandomLevel( unit = unique(studyDesign$transect) )
# Random Structure for Year (temporal data)
td <- data.frame(year=unique(M$Year))
row.names(td) <- unique(M$Year)
rL_year = HmscRandomLevel( sData = td )
# Random structure for sites
rL_site <-  HmscRandomLevel( unit = unique(studyDesign$site) )
# transect year
rL_ty <- HmscRandomLevel( unit = unique(studyDesign$ty) )





## formula for fixed effects
XFormula = ~ poly(shore.height, degree = 2, raw = TRUE)*year
  # poly(shore.height, degree = 2, raw = TRUE):year
  # anom.pine.sum.1 + anom.pine.win
  
  

## Traits -- will include the kelp.fucoid.turf type of functional grouping 
#            from the file "Algae_functional_groups.csv", but restricted to algae only

## Phylogeny -- may only be able to do this for red algae or brown algae separately


m <- Hmsc( Y = Y, 
           XData = XData, XFormula = XFormula,
           distr = "poisson",
           studyDesign = studyDesign, 
           ranLevels = list(site=rL_site, transect=rL ) )#,
                            # year=rL_year, 
                             # ty=rL_ty ) )  

# save the prepared model so it can be run in a separate script
save( m, file="R Code and Analysis/output from r/hmsc_specified.Rdata" )



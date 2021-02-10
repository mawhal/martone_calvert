# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script uses HMSC to model communities change over time and space

# Running HURDLE model
# see https://www.helsinki.fi/en/researchgroups/statistical-ecology/hmsc
  # Wednesday 4th November - R-demonstration 4 youtube video (39minutes)

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
  mutate( taxon = gsub(" ",".",taxon_lumped3) ) %>% 
  group_by( UID, Year, Site, Zone, taxon ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) 

# # average cover per transect
# dmean <- d.simple %>%
#   spread( taxon_lumped2, Abundance, fill=0 ) %>%
#   gather( taxon_lumped2, Abundance, -UID, -Year, -Site, -Zone ) %>%
#   group_by( Year, Site, Zone, taxon_lumped2 ) %>%
#   summarise( Abundance=mean(Abundance) )

d.comm.prep <- d.simple  # dmean

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.comm.prep %>%
  spread( taxon, Abundance, fill=0 )

# order community data by site and zone
d.comm <- d.comm %>%
  arrange( Site, factor(Zone,levels=c("LOW","MID","HIGH")) )
d.comm$Zone <- factor( d.comm$Zone,levels=c("LOW","MID","HIGH") )



# isolate the community, site, and sample data
comm.all <- d.comm[,-c(1:4)]
comm.all <- ceiling(comm.all)

# 
windows(12,3)
par( mar=c(3,4,0.5,0.5)+0.01, cex=0.7, las=1, cex=1.1 )
# boxplot( comm.all[rev(order(colSums(comm.all)))], pch=16, cex=0.3, axes=F )
# axis(2)
# axis(1, at=c(1,seq(10,300,by=10)) )
# all data
boxplot( comm.all[rev(order(apply(comm.all,2,function(z) length(z[z>0]))))], 
         pch=16, cex=0.3, axes=F, col = 'moccasin', outwex = 0.5, boxcol = 'black', outcol = scales::alpha("black",0.25) )
# # non-zero cover
# comm.present <- comm.all
# comm.present[ comm.present == 0] <- NA
# boxplot( comm.present[rev(order(apply(comm.all,2,function(z) length(z[z>0]))))],
#          pch=16, cex=0.3, axes=F, col = 'moccasin', log='y' )
axis(2)
mtext("cover (%)",2,line = 3,las=3)
axis(1, at=c(1,seq(10,300,by=10)) )
abline(v=48,col='red')
# abline(h=5,lty=1)
mtext("occupancy rank",1,line = 2)
box()
dev.off()

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
abundance <- colSums(comm)
which( colSums(comm) <=10 )
# comm <- comm[ ,-which( abundance <=10 )  ]
# comm <- comm.all
commpa <- comm
commpa[commpa>0] <- 1
sort(colSums(commpa))[1:20]
# reorder community matrix based on the most abundant to least (across all quadrats)
Y <- as.matrix( comm[, order(colSums(comm),decreasing = T) ] )
# convert Y to presence-absence matrix, and then to abundance, given presence (zeros become NA)
Ypa <- ifelse( Y>0,1,0 )
Yap <- ifelse( Y==0, NA, Y )
# log transform abundance given presence, so we can run a gaussian model
Ylap <- log( Yap )

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

# look at genus overlap
sort(unlist(lapply( strsplit( colnames(Y), split = "[.]"), function(z) z[1] )))


# define the variables to test from metadata and data
# merge community data with metadata
metacomm <- left_join( d.comm, muse )
write_csv( metacomm, "R Code and Analysis/output from r/community.csv")

# # read temperature data
# pine <- read_csv( "R Code and Analysis/output from r/PineIsland_summary.csv" )

# merge temperature and the rest of the metadata
# M <- left_join( metacomm, pine )
M <- metacomm
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
XData.raw <- M %>%
  ungroup() %>%
  select( year=Year, shore.height = Shore_height_cm )#, anom.pine.sum.1, anom.pine.win ) #%>%
  #mutate( year=factor(year, ordered=TRUE ) )
XData.raw <- as.data.frame( XData.raw )

# centered and scaled 
XData <- data.frame( XData.raw$year,XData.raw$shore.height, as.data.frame(poly(XData.raw$year,2)),
            as.data.frame(poly(XData.raw$shore.height,2)) )
names(XData) <- c("year","elev","year1","year2","elev1","elev2")

# set up random effects for hmsc
# make sure to use spatial data
# STUDY DESIGN
studyDesign <- data.frame( year = factor(M$Year),
                           quadrat = factor(M$UID),
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
# qudarat level (based on combination of site, transect, and quad)
rL_quad <- HmscRandomLevel( unit = unique(studyDesign$quadrat) )




## formula for fixed effects
XFormula = ~ poly(shore.height, degree = 2, raw = F)*poly(year, degree=5, raw=F)
XFormula = ~ year1 + year2 + elev1 + elev2 +
  elev1:year1 + elev1:year2 #+ elev2:year1  
  
  

## Traits -- could include the kelp.fucoid.turf type of functional grouping 
#            from the file "Algae_functional_groups.csv", but restricted to algae only
#           could also include generally occpuied zone (e.g., high, mid, low)

## Phylogeny -- may only be able to do this for red algae or brown algae separately


# For hurdle model, we need to run two model
# probit presence-absence model with ones and zeros
mprobit <- Hmsc( Y = Ypa, 
           XData = XData, XFormula = XFormula,
           distr = "probit",
           studyDesign = studyDesign, 
           ranLevels = list(site=rL_site, transect=rL) ) #, quadrat=rL_quad ) )#,
                            # year=rL_year, 
                             # ty=rL_ty ) ) 
# second, Gaussian model with log-abundances given presence
mnormal <- Hmsc( Y = Ylap, YScale = TRUE,
                 XData = XData, XFormula = XFormula,
                 distr = "normal",
                 studyDesign = studyDesign, 
                 ranLevels = list(site=rL_site, transect=rL) )

# combine these models into a list
models = list( mprobit, mnormal )
modelnames= c( "presence_absence", "abundance_COP" )

# save the prepared model so it can be run in a separate script
# save( mprobit, file="R Code and Analysis/output from r/hmsc_hurdle_probit_specified.Rdata" )
# save( mnormal, file="R Code and Analysis/output from r/hmsc_hurdle_normal_specified.Rdata" )
save(models, modelnames, file="R/output/hmsc_hurdle_specified.Rdata")


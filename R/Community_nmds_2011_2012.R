# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script identifies common/indicator taxa of the different elevation zones across sites 

# load libraries
library( tidyverse )
library( vegan )
#


## read data files
# all data that has been cleaned, taxon names corrected, and with lumping names and functional groups
ad <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_data_lump_function.csv" )
# all metadata
am <- read_csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )


## Data cleaning for Analysis -- consider moving part of this to another script
# use only 2011 and 2012 data, closest to establishment of transects
muse <- am[ am$Year %in% 2011:2012, ]

# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]
# Only use Mid-shore transects for now
# muse <- muse[ muse$Zone == "MID", ]
# muse <- muse[ muse$Site == "North Beach", ]
muse <- droplevels(muse)

# restrict rows of d to ones with UID's in the metadata
duse <- ad[ ad$UID %in% muse$UID, ]
dm <- left_join( duse, muse )

dm %>% 
  select(Site,Year,Quadrat) %>% 
  distinct() %>% 
  table()


# choose which to include
d <- dm %>% 
  filter( motile_sessile == "sessile" ) #%>% 
  # filter( non.alga.flag %in% c("Algae")  ) # restrict community analysis to algae only

# take a closer look
d %>% 
  group_by(taxon_lumped3) %>% 
  summarize(total=sum(Abundance)) %>% 
  arrange( total )

# average elevation per zone
dmeanelev <- d %>% 
  group_by( Year, Site, Zone ) %>%
  summarise( Elevation=mean(Shore_height_cm, na.rm=T) )
  
# calculate mean abundance per transect in each year
# first fix some taxa names
d$taxon_lumped[d$taxon_lumped=="Mytilus trossulus"] <- "Mytilus sp."
d$taxon_lumped2[d$taxon_lumped2=="Haliclona"] <- "Tunicata/Porifera"
d$taxon_lumped2[d$taxon_lumped2=="Halichondria"] <- "Tunicata/Porifera"

# remove rare, poorly resolved, or inconsistently counted taxa 
d <- d %>% 
  filter( !(taxon_lumped %in% c("articulated coralline","Unknown red blade",
                                "Colonial Diatoms", "Pleonosporium vancouverianum",
                                "Acrochaetium sp.", "Fauchea") ) )



## taxa not in final list
# remove "Chiharaea bodegensis", "Desmarestia herbacea", "Ectocarpus siliculosus"
# remove "Gloiocladia sp.", "Kornmannia leptoderma", "Petalonia fascia", "Tiffaniella snyderae", "Rhodymenia sp."
# "Phycodrys sp." = Polyneura latissima
d <- d %>% 
  filter( !(taxon_lumped3 %in% c("Chiharaea bodegensis", "Desmarestia herbacea", "Ectocarpus siliculosus",
                                 "Ulothrix-Urospora sp.") ) )
d <- d %>% 
  filter( !(taxon_lumped3 %in% c("Gloiocladia sp.", "Kornmannia leptoderma", "Petalonia fascia", "Tiffaniella snyderae", "Rhodymenia sp.") ) )
sort(unique(d$taxon_lumped3))

# add together taxa that are not unique to each quadrat
# this uses lumped taxon names, which will reduce the size of the dataset a bit
# restrict this to seaweeds and sessile invertebrates
d.simple <- d %>%
  mutate( taxon = gsub(" ",".",taxon_lumped3) ) %>% 
  group_by( UID, Year, Site, Zone, Quadrat, Meter.point, taxon=taxon_lumped3, funct_2021 ) %>%
  summarize( Abundance=sum(Abundance,na.rm=T)) 
d.simple %>%
  ungroup() %>% 
  select(Site,Year,Quadrat) %>% 
  distinct() %>% 
  table()

# pivot wider
dwide <-  d.simple %>%
  group_by( Year, Site, Zone, Quadrat, taxon ) %>%
  summarise( Abundance=mean(Abundance) ) %>% 
  pivot_wider( names_from = taxon, values_from = Abundance, values_fill=0 ) %>% 
  ungroup()

dmean <- d %>% 
  group_by( Year, Site, Zone, taxon = taxon_lumped3 ) %>%
  # filter( taxon %in% taxa_keep ) %>% 
  summarise( Abundance=mean(Abundance) )


# spread out
d.comm.mean <- dmean %>%
  pivot_wider( names_from = taxon, values_from = Abundance, values_fill=0 ) %>%
  ungroup() %>% 
  mutate( Zone = factor(Zone, levels=c("LOW","MID","HIGH")) ) %>% 
  arrange( Year, Site, Zone )


# remove UID column from community data
meta <- d.comm.mean[ ,1:3 ]
meta <- left_join(meta,dmeanelev)
comm <- as.matrix(d.comm.mean[,-c(1:3)])
dim(comm)
# interrogate the dataset
sort( colSums(comm), decreasing = T )


# NMDS
mds <- metaMDS(comm,
        distance = "bray",
        k = 3,
        maxit = 999, 
        trymax = 500,
        wascores = TRUE)

meta$Zone <- factor( meta$Zone, levels = c("LOW","MID","HIGH"))
meta$Site 
plot(mds)
ordihull(
  mds,
  meta$Zone,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("black", "black", "gray48" ),
  lty = c(1, 2, 1),
  lwd = 2.5
)
# strong separation of sites. Which taxa are most indicative of zones?


# Simper analysis
simp <- simper( comm, group = meta$Zone )
# grab top 5 taxa in each category
top <- lapply(simp, function(z) sort(z$average,decreasing = T)[1:10])
names(top) <- c("","","")

# names of the species of interest
topnames <- unique(names(do.call(c, top )))

Y <- scores(mds, display="species")
Y <- Y[rownames(Y) %in% topnames,]


# plot
par(mar=c(5,4,2,2)+0.1 )
plot(mds, type="n", ylim = c(-1,1),xlim = c(-1.5,1.5))
points( scores(mds)$sites[meta$Site=="Foggy Cove",1:2], display = "sites", 
        pch=21, bg="black" )
points( scores(mds)$sites[meta$Site=="Fifth Beach",1:2], display = "sites", 
        pch=21, bg="red" )
points( scores(mds)$sites[meta$Site=="North Beach",1:2], display = "sites", 
        pch=21, bg="blue" )
text( x = Y[,1], y = Y[,2], rownames(Y), cex = 0.5 )
ordihull(
  mds,
  meta$Zone,
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("black", "black", "gray48" ),
  lty = c(1, 2, 1),
  lwd = 2.5
)

### Martone seaweed data from Calvert Island, Central Coast BC
# This script reads and exports data into a community dataset (rows
# are samples and columns are taxa) for each survey year 
# 
# started by Matt Whalen on 8 sep 2021
#

library(tidyverse)

# all metadata
am <- read.csv( "data/R Code for Data Prep/Output from R/Martone_Hakai_metadata.csv", stringsAsFactors = TRUE )

## Data cleaning for Analysis -- consider moving part of this to another script
# remove 2011 data
muse <- am[ am$Year != "2011", ]
# remove Meay Channel
## NOTE THAT THIS ANALYSIS DOES NOT REQUIRE EQUAL SAMPLING OVER TIME OR SPACE
## a mjor exception to this is for convergence analysis, see trajectoryConvergence()
muse <- muse[ muse$Site != "Meay Channel", ]

# KEEP all quads, even if elevation is unknown
# # no NA values allowed, so we need to remove these from the dataset
# rem <- unique( which(is.na(muse$Shore_height_cm))  )
# # get rid of row 2 in all data structures
# muse   <- muse[-rem,]
muse <- droplevels(muse)


# read dataset generated in script "Community_rda.R"
d.simple <- read_csv("R/output/data_select_rda_HMSC.csv")
# Note that some taxa have already been removed. Please ask Whalen about these
# if you are interested

# spread Taxon column out into many columns filled with abundance/cover data
d.comm <- d.comm.prep %>%
  select( -funct_2021 ) %>%
  spread( taxon, Abundance, fill=0 )

# order community data by site and zone
d.comm <- d.comm %>%
  arrange( Site, factor(Zone,levels=c("LOW","MID","HIGH")) )
d.comm$Zone <- factor( d.comm$Zone,levels=c("LOW","MID","HIGH") )



# optional to deal with cover of 0.5
# comm.all <- ceiling(comm.all)


# define the variables to test from metadata and data
# merge community data with metadata
metacomm <- left_join( d.comm, muse )
write_csv( metacomm, "R/output/community.csv")

# get rid of extra columns
commselect <- d.comm[, -c(2,3,4,5)]

# isolate a data set for each year
commsplit <- split( commselect, f = metacomm$Year )

# write these to disk
sapply( names(commsplit), 
       function (x) write_csv( commsplit[[x]], file = paste0("R/output/community_", x, ".csv") )   )

x = "2012"
write_csv( as.data.frame(commsplit[[x]]), file = paste0("output/community_", x, ".csv") )





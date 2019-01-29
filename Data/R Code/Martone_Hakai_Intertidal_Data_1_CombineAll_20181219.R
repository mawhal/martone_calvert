###################################################
### Martone Lab Rocky Intertidal Community Data ###
###     All Data collected on Calvert Island    ###
###################################################

### The purpose of this script is to join all
### available data at it is organized in its 
### rawest form

###################################################
### Data are percent cover of seaweeds from 
###   quadrats along transects in low, mid,
###   and high tide zones
###################################################
###
### code created by Matt Whalen, 3 April 2018
### updated 19 December 2018


# load libraries
library(tidyverse)
library(readxl) 
library(reshape2)




# NOTES
# Fix typos in Abundance - look for commas


########################################################
## need to find a way to read all excel files at once, 
##    so that data entered in the future can easily be incoporated into the dataset
## readxl package (see https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames/12948450)
## see also http://readxl.tidyverse.org/articles/articles/readxl-workflows.html


## function to read in data using readxl 
read_excel_all <- function( data ){
  path <- paste0("../Excel Files/All years/edited/",data)
  sheets <- path %>% 
    excel_sheets() %>% 
    set_names() %>% 
    map( read_xlsx, path = path, col_names=FALSE  )
  return( sheets )
}

ptm2012 <- read_excel_all( "2012 Hakai_edited_MAW.xlsx" )
ptm2018 <- read_excel_all( "2018_Hakai_edited_MAW.xlsx" )



#### ISOLATE DATA AND METADATA (e.g. location info)
test <- ptm2018
# we now have a list, each element of which is a sheet from the xlsx file
# str(ptm2012)
names(test)
# Goal: create a single data.frame with rows being data from individual plots
# all sheets have same number of columns, but different number of rows
lapply( ptm2018, dim )

#
## METADATA
# metadata is sheet name and first X rows of data.frame
# number of rows by year: 2011 (5), 2012 (5), 2013 (8), 2014 (12), 2015 (12), 2016 (12), 2017 (12)

extractYEAR <- function( data, sheetnames, header=13 ){
  metarows <- t( data[1:header,] )                #transpose to flip columns to rows
    meta <- as_tibble( metarows[-1,],    # generate data.frame without column headings
                           stringsAsFactors = FALSE, row.names = NULL ) 
    colnames(meta) <- metarows[1,]           # add column headings
    # meta$SiteHeightYear <- sheetname # add sheet name as a column
  # NOTE: date is in excel format (will deal with this later once we have all metadata together)
  # NOTE: Meter Point column should have whole numbers, but does not always...need to round
  #
  ## DATA
  # data is everything else after the rows containing metadata
  draw <- data[header+1:nrow(data),]
    colnames(draw) <- c("Taxon",1:10)
    # melt the dataframe
    dmelt <- melt( draw, id.vars = 1, 
                   variable.name = "Quadrat", value.name = "Abundance") 
    # NOTE Abundance can mean percent cover or count or presence
    # get rid of NA's and 0's
    d <- na.exclude(dmelt)
    d <- d[d$Abundance>0,]
    # d$SiteHeightYear <- sheetname # add sheet name as a column, this will help us merge later
    d <- as_tibble(d)
  return( list(meta=meta,data=d) )
}

# extractYEAR( test )

# lapply the extract function to all sheets in the xlsx file
testsep  <- lapply( test, extractYEAR )


# a series of rbinds will group the data together into a meta and a data data.frame
testbind <- do.call( rbind, testsep )
# define names for the sheets
sheetnames  <- names(testsep)
# add sheet names as columns in each list element
testbindnames <- mapply( cbind, testbind, "SiteHeightYear"=sheetnames, SIMPLIFY = TRUE)

# separate metadata and data parts of the list
meta <- as_tibble(do.call( rbind, testbindnames[1:(length(testbindnames)/2)]))
d    <- as_tibble(do.call( rbind, testbindnames[(length(testbindnames)/2)+1:length(testbindnames)]))

list( meta=meta, data=d )


######## Wrap everything into a single function
#
# INTERATE OVER ALL SHEETS (wrap into a function, then lapply across sheets)
# all metadata, then all data, then merge the two together?

# lapply over all years of data
# list files
allyears <- list.files( pattern=".xlsx", path="../Excel Files/All years/edited/" )

# single function to take each file and run the extractYEAR function on every sheet,
# given a number of header rows
# number of rows by year: 2011 (10), 2012 (10), 2013 (8), 2014 (12), 2015 (12), 2016 (12), 2017 (12), 2018 (13)
header <- c( 10,8,12,12,12,12,13 )
extractALL <- function( files, header=header ){
  # read all of the files
  data = lapply( allyears, read_excel_all )
  # str( data, max.level=2 )
  # run the extractYEAR function for all sheets in all files
  metalist = list()
  dlist    = list()
  for( i in 1:length(data) ){
    # run extractYEAR on each data piece
    temp     = lapply( data[[i]], extractYEAR, header=header[i] )
    # a series of rbinds will group the data together into a meta and a data data.frame
    tempbind = do.call( rbind, temp )
    # define names for the sheets
    sheetnames  = names(temp)
    # add sheet names as columns in each list element
    tempbindnames = mapply( cbind, tempbind, "SiteHeightYear"=sheetnames, SIMPLIFY = TRUE)
    # separate metadata and data parts of the list
    meta = as_tibble(do.call( rbind, tempbindnames[1:(length(tempbindnames)/2)]))
    d    = as_tibble(do.call( rbind, tempbindnames[(length(tempbindnames)/2)+1:length(tempbindnames)]))
    metalist[[i]] = meta
    dlist[[i]]    = d
  }
  return( list(metalist,dlist) )
}


all.list <- extractALL( allyears, header )

# look at list structure
str( all.list, max.level = 2 )
# first element contains a list of all metadata. These have different numbers of columns as data recording changed over time
# second element contains cover/abundance/presence data. These do contain the same number of columns...can rbind

### rbind data
all.data <- do.call( rbind, all.list[[2]] )
all.data

### combine metadata
# extract metadata from list
meta.list <- all.list[[1]]
# use Reduce to successively join all of the metadata
all.meta <- Reduce( full_join, meta.list ) # some warnings are okay here
# format dates
all.meta$Date <- as.Date( as.numeric(all.meta$Date), origin = "1899-12-30")
# correct all spelling of Meay Channel
all.meta$SiteHeightYear <- gsub("Maey", "Meay", all.meta$SiteHeightYear)
# expand SiteHeightYear into three separate columns
splits <- strsplit( all.meta$SiteHeightYear, split=" " ) 
all.meta$Site   <- unlist( lapply( splits, function(z) paste( z[1], z[2], sep=" " ) ) )
all.meta$Zone <- unlist( lapply( splits, function(z) z[3] ) )
all.meta$Year   <- unlist( lapply( splits, function(z) z[4] ) )
# currently meter point has lots of different values
sort(unique( all.meta$`Meter point` ))
# get rid of spaces and comments
all.meta$`Meter point`[ all.meta$`Meter point` == "   30" ] <- 30
all.meta$`Meter point`[ all.meta$`Meter point` == "   45" ] <- 45
all.meta$`Meter point`[ all.meta$`Meter point` == "11 (instead of 12)" ] <- 11
all.meta$`Meter point`[ all.meta$`Meter point` == "16 (instead of 15)" ] <- 16
# make Meter point numeric and round it
all.meta$`Meter point` <- round( as.numeric( all.meta$`Meter point`) )



######################
### combine Elevation data with metadata
## NOTE: at the time of writing this script, elevation data were not available for the Maey Channel Site

# some surveys listed erroneous numbers for meters along transect
# fifth mid 2012, quad 10 is 26 meters, quad 9 is 25
# west low 2014, quad 9 is 35
# north low 2017, quad 1 is 1


# read elevation data
elev <- read_xlsx( "../Shore Heights Elevation/Elevation_transects2014_2.xlsx" )
# rename entries in Site and Zone
elev$Site <- plyr::revalue( elev$Site, c("North"="North Beach", "Fifth"="Fifth Beach", "West"="West Beach") )
# make zone designations all upper case
elev <- plyr::mutate( elev, Zone=toupper(Zone) )
# select and rename columns
elev <- elev %>%
  select( Site, Zone, 'Meter point'=Transect_num, Shore_height_cm )

all.meta <- left_join( all.meta, elev )

# write missing data to file
write.csv( all.meta[ is.na(all.meta$Shore_height_cm) & all.meta$Site!="Meay Channel", ], 
           "Output from R/Martone_Hakai_missingElevation.csv", row.names=F  )
all.meta[ is.na(all.meta$Shore_height_cm) & all.meta$Site!="Meay Channel", ]
# missing ones (n=3) are either at zero meters along the transect (n=1) 
    # or at distances along the transect that are greater than have been surveyed for elevation (n=2)


# Data are now essentially combined. There are two tables (tibbles) that contain 1) quadrat information and 2) data
# Note, data still needs cleaning before analysis
    # for example. Column abundance (which is usually % cover) has lots of character values instead of numbers
# In addition, taxa might change over the years, and they need consistent labels
    # in a different script, use the Master Species Codes spreadsheet to adjust names
# SiteHeightYear (the column used to merge data and metadata) also has some inconsistency over time
  # Matt manually changed this for the excel files

# Save these and consider this script done...for now
dim( all.data )
dim( all.meta )

write.csv( all.data, "Output from R/Martone_Hakai_data.csv", row.names=F )
write.csv( all.meta, "Output from R/Martone_Hakai_metadata.csv", row.names=F )
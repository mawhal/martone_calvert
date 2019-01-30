###################################################
### Martone Lab Rocky Intertidal Community Data ###
###     All Data collected on Calvert Island    ###
###################################################
# 
# by Matt Whalen
# updated 29 January 2019

# This script cleans up the data and deals with mistakes and taxonomic synonymy

# set options
options( stringsAsFactors = FALSE )

# load libraries
library(tidyverse)


## read data files
# all data
ad <- read.csv( "Output from R/Martone_Hakai_data.csv" )
# all metadata
am <- read.csv("Output from R/Martone_Hakai_metadata.csv" )


## Deal with trace cover and other oddities
# # replace all commas with periods for Abundance
# ad$Abundance <- as.numeric( gsub( ",", "[.]", ad$Abundance ) )
# change "present" to 0,5
ad$Abundance <- gsub( "present", "0.5", ad$Abundance )
# change 1 of Fucus to "trace"
ad$Abundance <- gsub( "1 on Fucus", "trace", ad$Abundance ) 
# change trace to 0.5% cover
ad$Abundance <- gsub( "t.*", "0.5", ad$Abundance, ignore.case = TRUE ) 
# change "0,5" to "0.5"
ad$Abundance <- gsub( "0,5", "0.5", ad$Abundance, ignore.case = TRUE ) 

# accept only the first thing if separated by certain characters
asplit <- strsplit( ad$Abundance, split =c("/|;"))
ad$Abundance <- unlist( lapply( asplit, function(z) z[1] ))
# LOB = lots of babies
ad$Abundance <- gsub( "LOB","100",ad$Abundance )
# 10s and 100s
ad$Abundance <- gsub( "10s","20",ad$Abundance )
ad$Abundance <- gsub( "100s","200",ad$Abundance )
# spaces and parentheses
ad$Abundance <- gsub( "[ (].*","",ad$Abundance )

# look at uniue entries and their conversions
uni <- sort(unique(ad$Abundance))
data.frame(character=uni, number=as.numeric(uni))
ad$Abundance <- as.numeric( ad$Abundance )



## Taxon names
# get rid of things after and including parentheses
ad$Taxon <- gsub( " [(].*", "", ad$Taxon )
# get rid of percent signs
ad$Taxon <- gsub( "%", "", ad$Taxon )
# get rid of plus signs
ad$Taxon <- gsub( "[+].*" , "", ad$Taxon )
# get rid of quotes
ad$Taxon <- gsub( '".*',"", ad$Taxon )
# get rid of double dash
ad$Taxon <- gsub( "--"," ", ad$Taxon )
# get rid of single dash
ad$Taxon <- gsub( " - "," ", ad$Taxon )

# get rid of with or on
ad$Taxon <- gsub( "with.*", "", ad$Taxon)

# for the sake of reducing unique entries, capitalize first letter of all Taxon names
ad$Taxon <- paste0( toupper(substr(ad$Taxon, 1, 1)), substr(ad$Taxon, 2, nchar(ad$Taxon)) )



# taxa specific things, e.g. fixing mispellings
ad$Taxon <- gsub( "Acroch.*", "Acrochaetium sp.", ad$Taxon )
ad$Taxon <- gsub( "Amphipod.*", "Amphipoda", ad$Taxon )
ad$Taxon <- gsub( ".*diatoms.*", "diatoms", ad$Taxon )
ad$Taxon <- gsub( "Kathrina", "Katharina", ad$Taxon )
ad$Taxon <- gsub( "Katerina", "Katharina", ad$Taxon )
ad$Taxon <- gsub( "Pagarus", "Pagurus", ad$Taxon )
ad$Taxon <- gsub( "Panicea", "Halichondria panicea", ad$Taxon )
ad$Taxon <- gsub( "Lottornia", "Littorina", ad$Taxon )
ad$Taxon <- gsub( "Lottorine snails", "Littorina", ad$Taxon )
ad$Taxon <- gsub( "subsp ", "subsp. ", ad$Taxon )
ad$Taxon <- gsub( "Anthopleura ele.*","Anthopleura elegantissima", ad$Taxon )
ad$Taxon <- gsub( "Anthopleura x.*","Anthopleura xanthogrammica", ad$Taxon )
ad$Taxon <- gsub( "Anthopleura z.*","Anthopleura xanthogrammica", ad$Taxon )
ad$Taxon <- gsub( "bare.*","bare rock", ad$Taxon, ignore.case = TRUE )
ad$Taxon <- gsub( "Barnacles, goose-necked","Pollicipes", ad$Taxon )
# genus "M." appears to be Mazaella
ad$Taxon <- gsub( "M[.]", "Mazzaella ", ad$Taxon )
# Species "M.lat is probably Mazzaella latissimus
ad$Taxon <- gsub( "Mazzaella lat.*", "Mazzaella latissimus", ad$Taxon )
ad$Taxon <- gsub( "Lophopan", "Lophopanopeus bellus", ad$Taxon )
ad$Taxon <- gsub( "Pagurs", "Pagurus", ad$Taxon )
ad$Taxon <- gsub( "Holiclona", "Haliclona", ad$Taxon )
# plural to singular
ad$Taxon <- gsub( "Crabs", "Crab", ad$Taxon )
ad$Taxon <- gsub( "Snails", "Snail", ad$Taxon )
ad$Taxon <- gsub( "Mussels", "Mussel", ad$Taxon )
ad$Taxon <- gsub( "Limpets", "Limpet", ad$Taxon )
ad$Taxon <- gsub( "Kelp Crab", "Kelp crab", ad$Taxon )
ad$Taxon <- gsub( "Isopods", "Isopoda", ad$Taxon )
ad$Taxon <- gsub( "Hydroids", "Hydroid", ad$Taxon )

# for numbered species, make a rule about how it should look
ad$Taxon <- gsub( "sp([0-9]+)", "sp.\\1", ad$Taxon )

# trim all the white space
ad$Taxon  <- trimws( ad$Taxon )



# "taxa" to remove
ad <- ad[ !ad$Taxon %in% c( "Animal notes", "Animals" ), ]
ad <- ad[ ad$Taxon != c( "Habitat notes" ), ]

# consider using tidyverse filter to make some of the above code a little prettier


## Use corrected species names to replace taxon names
corrected_taxa <- read.csv( "../taxa/CorrectedTaxonList.csv" )
ad.corrected <- left_join( ad, corrected_taxa, by=c("Taxon"="taxon") )

# fix the ones that weren't covered
unique(ad.corrected$Taxon[ is.na(ad.corrected$taxon_corrected) ])
ad.corrected$taxon_corrected[ ad.corrected$Taxon == "Lophopanopeus bellus" ] <- "Lophopanopeus bellus"
ad.corrected$taxon_corrected[ ad.corrected$Taxon == "bare rock" ] <- "bare rock"
unique(ad.corrected$Taxon[ is.na(ad.corrected$taxon_corrected) ])

sort(unique(ad$Taxon))
sort(unique(ad.corrected$taxon_corrected))

# recreat ad with corrected names, then sum across corrected taxa
ad <- ad.corrected %>%
  select( UID, Taxon=taxon_corrected, Abundance ) %>%
  group_by( UID, Taxon ) %>%
  summarise( Abundance=sum(Abundance) )

# Cobble. need to move this to metadata!
ad[ad$Taxon=="Cobble",]
ad[ad$Taxon=="BARE ROCK / SUBSTRATE (%)",]


# save the data to disk, overwriting the previous datafile
write.csv( ad, "Output from R/Martone_Hakai_data.csv", row.names = FALSE )

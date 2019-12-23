# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen

# This script generates a table of occurrences across time

# load libraries
library( tidyverse )


# load the data
metacomm <- read_csv( "R Code and Analysis/output from r/community.csv" )

# replace space with period
names(metacomm) <- gsub( " ", ".", names(metacomm) )

# use actual data
comm_select <- metacomm %>% 
  gather( key = taxon, value = N, -UID ) %>%
  filter( N>0 )
comm_select <- comm_select %>% 
  separate( UID,c("beach","b","zone","year","quad") ) %>% 
  unite( site, beach, b )
# # remove rare taxa
# comm_rare <- comm_select %>% 
#   group_by( taxon ) %>% 
#   summarize( N=sum(N) ) %>% 
#   filter( N<=10 )
# # merge
# comm_final <- comm_select %>% 
#   filter( !(taxon %in% comm_rare$taxon) )


inverts <- c("Barnacles","Mytilus.sp.","Anemone","Bryozoan","Tunicata/Porifera","Pollicipes.polymerus","Tube.worms","Hydroid" )
comm_final$alga <- "alga"
comm_final$alga[ comm_final$taxon %in% inverts ] <- "invert"



# occurrences across years
commtab <- with(comm_select, table(taxon, year))
commtabdf <- as.data.frame.matrix( commtab )
occurrence <- data.frame( taxon=rownames(commtabdf),commtabdf)
# commtab <- commtab[ order(rowSums(commtab),decreasing = T), ]
write_csv( occurrence, "R Code and Analysis/output from r/occurrence_table.csv")

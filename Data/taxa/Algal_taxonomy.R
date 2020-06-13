### grab taxonomy for algae

# libraries
library(tidyverse)
library(taxize)

# read the data
algae <- read_csv("Data/taxa/CorrectedTaxonList.csv") %>% 
  filter( kingdom=="Algae")
# take first name
first <- unlist(lapply( strsplit(algae$taxon_corrected,split = " "), 
                        function(z) z[1] ))
names <- sort(unique(first))

tax <- classification( names, db = 'ncbi' )
taxbind <- rbind(tax)
taxspread <- taxbind %>% 
  select(-id) %>% 
  filter( rank != "no rank" ) %>% 
  group_by(query) %>% 
  spread( rank, name )
# fill in brown algal phylum
taxspread$phylum[ taxspread$class == "Phaeophyceae" ] <-  "Ochrophyta"


# write to disk
write_csv( taxspread, "Data/taxa/algal_taxonomy.csv")

# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
#
# started 23 June 2022
#
# this script identifies the dominant taxa in each transect in each year
# to help inform questions about changes in zonation patterns over time

# packages
library(tidyverse)
library(vegan)

# read community data
comm.all <- read_csv("R/output/community_all.csv" )

# isolate community data and metadata
comm <- comm.all %>% select(Acrosiphonia:Unknown.red.blade)
metacomm <- comm.all[,c(1:4)]

## approaches 
# compare rank-abundance curves?
# simple top X number based on relative abundance

## calculate relative abundance
rel.abun <- decostand(comm, method = "total")
# for each row get the column order for relative abundance
domby <- by( as.matrix(rel.abun), factor(1:nrow(rel.abun)), order, decreasing = T )
domorder <- do.call(rbind, domby)
dom <- matrix( colnames(comm)[domorder], ncol = ncol(domorder) )

rel.comm <- bind_cols( metacomm, rel.abun )

# another option -  h.t. <https://stackoverflow.com/questions/62658868/how-to-find-dominant-species-per-plot>
dom.taxa <- rel.comm %>%
  select( UID, Year, Site, Zone, Acrosiphonia:Unknown.red.blade ) %>% 
  pivot_longer(cols = Acrosiphonia:Unknown.red.blade) %>%
  group_by(Year, Site, Zone, name) %>%
  summarize( value=mean(value) ) %>%
  # group_by(UID) %>% 
  arrange(desc(value)) %>%
  dplyr::mutate(cum_value = cumsum(value)) %>%
  slice(1:min(which(cum_value >= 0.5))) %>%
  dplyr::summarise(Species = paste(name, collapse = " "))

write_csv( dom.taxa, "R/output/dominant_taxa_transect_50percent.csv")

df <- structure(list(Plot = c("Plot_1", "Plot_2", "Plot_3"), Species_A = c(50L, 
                                                                           20L, 85L), Species_B = c(35L, 30L, 5L), Species_C = c(10L, 40L, 
                                                                                                                                 15L), Species_D = c(5L, 10L, 0L)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                        -3L))
df %>%
  pivot_longer(cols = starts_with("Species")) %>%
  group_by(Plot) %>%
  arrange(desc(value)) %>%
  dplyr::mutate(cum_value = cumsum(value)) %>%
  slice(1:min(which(cum_value >= 80))) %>%
  dplyr::summarise(Species = paste(name, collapse = " "))

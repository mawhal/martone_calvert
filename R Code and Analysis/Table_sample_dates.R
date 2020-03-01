# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen (adapted in part from Sam Starko's diversity script)
# created 26 Feb 2019


# load libraries
library(tidyverse)
library(vegan)
library(psych)
library(plotrix)


## read data files
# all metadata
am <- read_csv( "Data/R code for Data Prep/Output from R/Martone_Hakai_metadata.csv" )


# sampling dates in different years
View( am %>% 
  # mutate( month=)
  group_by( Date ) %>% 
  summarize(  samples=length(unique(UID)) ) )

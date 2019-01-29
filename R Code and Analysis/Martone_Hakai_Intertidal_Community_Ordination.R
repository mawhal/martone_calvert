# Martone Hakai Rocky Shore Seaweed Surveys
# 
# by Matt Whalen
# updated 20 December 2018

# This script runs ordination procedures on data

# set options
options(stringsAsFactors = FALSE)

# load libraries
library(tidyverse)


## read data files
# all data
ad <- read.csv( "../Data/Excel Files/All years raw/Output from R/Martone_Hakai_data.csv")
# all metadata
am <- read.csv("../Data/Excel Files/All years raw/Output from R/Martone_Hakai_metadata.csv" )


#